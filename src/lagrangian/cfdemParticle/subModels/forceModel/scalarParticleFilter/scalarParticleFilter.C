/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "scalarParticleFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(scalarParticleFilter, 0);

addToRunTimeSelectionTable
(
    forceModel,
    scalarParticleFilter,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
scalarParticleFilter::scalarParticleFilter
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    scalarTransportProperties_                  //this is clumsy, but effective
    (
        IOobject
        (
            "scalarTransportProperties",
            sm.mesh().time().constant(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    transportModelName_(propsDict_.lookup("scalarTransportModelName")),
    transportModelDict_(scalarTransportProperties_.subDict(transportModelName_+"Props")),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),             //common names/data
    velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
    partDat_(NULL),
    partDatFlux_(NULL),
    partDatTmpExpl_(NULL),
    partDatTmpImpl_(NULL),
    validPartFlux_(false),
    useModel2_(false),
    useDepositionRate_(true),
    useMomentumTransferRate_(false),
    eulerianFieldNames_( transportModelDict_.lookup("eulerianFields")), 
    isDropletCloud_(propsDict_.lookup("isDropletCloud")),
    partSpeciesNames_(propsDict_.lookup("partSpeciesNames")),
    partSpeciesFluxNames_(propsDict_.lookup("partSpeciesFluxNames")),
    partHeatFluxPositionInRegister_(-1),
    maxSource_(1e30),
    dropletDensity_(1000), //rhoLiq_("rhoLiq", dimMass/dimLength/dimLength/dimLength, 1000),
    dropletDiameter_(20e-6) //spray droplet diamter
{
    //Allocate and push back the pointer for transfer to external Code (e.g., LIGGGHTS)
    allocateMyArrays();
    if (propsDict_.found("partHeatFluxName"))
        particleCloud_.registerNamesFieldsUserCFDEMToExt(propsDict_.lookup("partHeatFluxName"),partHeatFluxPositionInRegister_);
 
    particleCloud_.checkAndregisterNamesFieldsUserCFDEMToExt(partSpeciesFluxNames_,partSpeciesFluxPositionInRegister_);

    //read other properties
    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }

    if (propsDict_.found("useModel1"))
    {
        useModel2_=readBool(propsDict_.lookup ("useModel1"));
        Info << "setting for useModel1: " << useModel2_ << endl;
    }
    if(useModel2_)
        LambdaParticle=&scalarParticleFilter::LambdaParticleModel2;
    else
        LambdaParticle=&scalarParticleFilter::LambdaParticleModel1;

    //Main switches
    if (propsDict_.found("useDepositionRate"))
    {
        useDepositionRate_=readBool(propsDict_.lookup ("useDepositionRate"));
        Info << "setting for useDepositionRate: " << useDepositionRate_ << endl;
    }
    if (propsDict_.found("useMomentumTransferRate"))
    {
        useMomentumTransferRate_=readBool(propsDict_.lookup ("useMomentumTransferRate"));
        Info << "setting for useMomentumTransferRate: " << useMomentumTransferRate_ << endl;
    }


    if(partSpeciesFluxNames_.size()>0)
    {
        Info << "scalarParticleFilter::partSpeciesFluxName(s): " << partSpeciesFluxNames_ << endl;
        validPartFlux_ = true;
    }

    if(!validPartFlux_)
        FatalError <<"You must set a valid heat flux name, or a valid transfer coefficient and fluid name \n" 
                   << abort(FatalError);

    dropletDensity_  =   scalar(readScalar(propsDict_.lookup("dropletDensity")));
    dropletDiameter_ =   scalar(readScalar(propsDict_.lookup("dropletDiameter")));

    particleCloud_.checkCG(true);

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch
    forceSubM(0).setSwitchesList(9,true); // activate verboseToDisk switch

    //set default switches (hard-coded default = false)
    //forceSubM(0).setSwitches(XXX,true);

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    // setup required communication
/*
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).setupCommunication();
*/
    Info << "scalarParticleFilter found the following speciesFieldNames: " << eulerianFieldNames_ << endl;

    // suppress particle probe
    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, typeName+".logDat");  
        particleCloud_.probeM().vectorFields_.append("Urel");               //first entry must the be the vector to probe
        particleCloud_.probeM().scalarFields_.append("Rep");                //other are debug
        particleCloud_.probeM().scalarFields_.append("alpha");              //other are debug
        particleCloud_.probeM().scalarFields_.append("voidfraction");       //other are debug
        particleCloud_.probeM().scalarFields_.append("LpRate");
        particleCloud_.probeM().writeHeader();
    }
    
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).constructorCalls(typeName);

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

scalarParticleFilter::~scalarParticleFilter()
{
    particleCloud_.dataExchangeM().destroy(partDat_,1);
    particleCloud_.dataExchangeM().destroy(partDatTmpExpl_,1);
    particleCloud_.dataExchangeM().destroy(partDatTmpImpl_,1);

    //external particle data (e.g., fluxes) will be destroyed in cloud.
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void scalarParticleFilter::allocateMyArrays() const
{
    if(particleCloud_.numberOfParticlesChanged())
    {
        //Only allocate temporary Lagrangian arrays, correct length will be provided by ExternalCode
        double initVal = 0.0;
        particleCloud_.dataExchangeM().allocateArray(partDat_,initVal,1);  
        particleCloud_.dataExchangeM().allocateArray(partDatTmpExpl_,initVal,1);  
        particleCloud_.dataExchangeM().allocateArray(partDatTmpImpl_,initVal,1); 
    }

    //external particle data (e.g., fluxes) will be allocated in cloud, 
    //just need to set pointers before access

}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void scalarParticleFilter::setForce() const
{
   if(!useMomentumTransferRate_) return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void scalarParticleFilter::manipulateScalarField(volScalarField& explicitEulerSource,
                                                 volScalarField& implicitEulerSource, //TODO: currently zero, split up similar to scalarGeneralExchange
                                                 int speciesID) const
{
    // reset Scalar field (== means hard reset)
    explicitEulerSource == dimensionedScalar("zero", explicitEulerSource.dimensions(), 0.);
    implicitEulerSource == dimensionedScalar("zero", implicitEulerSource.dimensions(), 0.);

    //Set the names of the exchange fields
    word    fieldName;
    word    partDatName;
    word    partFluxName;

    if(speciesID<0) //have temperature
    {
        if(partHeatFluxPositionInRegister_>-1)
            partDatFlux_           = particleCloud_.particleDatFieldsUserCFDEMToExt[partHeatFluxPositionInRegister_];
        else
            return;

        FatalError << "scalarParticleFilter does not implement fluxes for heat equation yet. You must NOT specify 'partHeatFluxName' in the sub-dictionary for the scalarParticleFiler." << abort(FatalError);
        if( forceSubM(0).verbose() )
        {
            Info <<"scalarParticleFilter will push particle heat fluxes (register:" << partHeatFluxPositionInRegister_ 
                 <<  "). \n" << endl;
        }
    }
    else
    {
        if(!isDropletCloud_[speciesID]) return;     //only perform on droplet cloud
        fieldName          = eulerianFieldNames_[speciesID];
        partDatName        = partSpeciesNames_[speciesID];
        partFluxName       = partSpeciesFluxNames_[speciesID]; 

        int positionInRegister = partSpeciesFluxPositionInRegister_[speciesID];
        partDatFlux_           = particleCloud_.particleDatFieldsUserCFDEMToExt[positionInRegister];

        if( forceSubM(0).verbose() )
        {
            Info <<"scalarParticleFilter manipulates eulerian field '"   << fieldName << "' now." << endl;
            Info <<"scalarParticleFilter will push particle fluxes to '" << partFluxName 
                 << "' (register:" << positionInRegister <<  "). \n" << endl;
        }
    }


    //==============================
    // get references to the current field that needs to be exchanged
    const volScalarField& voidfraction_(particleCloud_.mesh().lookupObject<volScalarField> (voidfractionFieldName_));    // ref to voidfraction field
    const volVectorField& U_(particleCloud_.mesh().lookupObject<volVectorField> (velFieldName_));
    const volScalarField& fluidScalarField_(particleCloud_.mesh().lookupObject<volScalarField> (fieldName));            // ref to scalar field
    //==============================

    // realloc the arrays
    allocateMyArrays();

    // get particle data
    if(partDatName!="none")
        particleCloud_.dataExchangeM().getData(partDatName,"scalar-atom", partDat_);

    const volScalarField& nufField  = forceSubM(0).nuField();
    const volScalarField& rhofField = forceSubM(0).rhoField();

      

    // calc La based heat flux
    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    scalar fluidValue(0);
    label  cellI=0;
    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar dParticle(0);
    scalar dparcel(0);
    scalar numberParticlesInParcel(1);
    scalar nuf(0);
    scalar rhof(0);
    scalar magUr(0);
    scalar As(0);
    scalar Rep(0); //mean slip Reynolds number
    scalar Stokes(0);
    scalar filterCoeff(0); //[1/m]

    #include "resetVoidfractionInterpolator.H"
    #include "resetUInterpolator.H"
    #include "resetFluidScalarFieldInterpolator.H"
    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if (cellI > -1) // particle Found
        {
            if(forceSubM(0).interpolation())
            {
                position = particleCloud_.position(index);
                voidfraction = voidfractionInterpolator_().interpolate(position,cellI);
                Ufluid = UInterpolator_().interpolate(position,cellI);
                fluidValue = fluidScalarFieldInterpolator_().interpolate(position,cellI);
            }else
            {
				voidfraction = voidfraction_[cellI];
                Ufluid = U_[cellI];
                fluidValue      = fluidScalarField_[cellI];
            }

            // calc relative velocity
            Us = particleCloud_.velocity(index);
            Ur = Ufluid-Us;
            magUr = mag(Ur);
            dParticle = 2.0*particleCloud_.radius(index);   //this is the true particle diameter!
            dparcel = dParticle;
            forceSubM(0).scaleDia(dParticle,index); //caution: this fct will scale ds!
            numberParticlesInParcel    = dparcel/dParticle;
            numberParticlesInParcel   *= numberParticlesInParcel*numberParticlesInParcel;
            
            As      = dParticle*dParticle*0.7853981634;                 //0.7853981634 = pi/4; projected surface area of the particle
            nuf     = nufField[cellI];                           
            rhof    = rhofField[cellI];
            Rep     = dParticle
                      * magUr
                      / nuf
                      * voidfraction;                 //TODO:  voidraction was missing here. please check again

            Stokes  = stokesNumber(magUr, dParticle, nuf*rhof);

            scalar lambdaParticle  = (this->*LambdaParticle)(Rep, Stokes, voidfraction);
            scalar alpha           = lambdaParticle  // single-collector efficiency *
                                   * fluidValue* voidfraction ;                         // particle concentration 

            scalar tmpPartFlux     = alpha* As * magUr/dropletDensity_;    //volumetric flux PER PARTICLE
            partDatFlux_[index][0]+= tmpPartFlux           //MUST ADD total source for ALL particles in parcel
                                   * numberParticlesInParcel;
            scalar tmpFlux = tmpPartFlux                    //total source for ALL particles in parcel
                            * numberParticlesInParcel;

            scalar areaTimesTransferCoefficient= -lambdaParticle 
                                                 * voidfraction  
                                                 * As
					     * magUr;

            forceSubM(0).explicitCorrScalar( partDatTmpImpl_[index][0], 
                                             partDatTmpExpl_[index][0],
                                             areaTimesTransferCoefficient,
                                             fluidValue,
                                             fluidScalarField_[cellI],
                                             partDat_[index][0],
                                             forceSubM(0).verbose()
                                           );

            if( forceSubM(0).verbose() )
            {


                filterCoeff  =  (this->*LambdaParticle)(Rep, Stokes ,voidfraction)  
                             *  1.5 
                             * (1.0-voidfraction) / dParticle;

                Pout << "index    = "   << index << endl;
                Pout << "partFlux = "   << tmpPartFlux << endl;
                Pout << "magUr = "      << magUr << endl;
                Pout << "As = "         << As << endl;
                Pout << "r = "          << particleCloud_.radius(index) << endl;
                Pout << "dParticle = "  << dParticle << endl;
                Pout << "nuf = "        << nuf << endl;
                Pout << "Rep = "        << Rep << endl;
            //  Pout << "lambdaP/Shp = " << (this->*LambdaParticle)(1.23456789,voidfraction) << endl;
                Pout << "lambdaP = "    << filterCoeff << endl;
                Pout << "voidfraction = " << voidfraction << endl;
                if(partDatName!="none")
                    Pout << "partDat_[index][0] = " << partDat_[index][0] << endl  ;
                Pout << "fluidValue = " << fluidValue << endl  ;
            }

            //Set value fields and write the probe; WARNING: MUST include "setupProbeModel.H" before adding data to the probe model!
            if(probeIt_)
            {
                #include "setupProbeModelfields.H"
                vValues.setSize(vValues.size()+1, Ur);
                sValues.setSize(sValues.size()+1, Rep);
                sValues.setSize(sValues.size()+1, alpha);
                sValues.setSize(sValues.size()+1, voidfraction);
                sValues.setSize(sValues.size()+1, tmpFlux);
                particleCloud_.probeM().writeProbe(index, sValues, vValues);
            }
        }
    }

    particleCloud_.averagingM().setScalarSum
    (
        explicitEulerSource,
        partDatTmpExpl_,
        particleCloud_.particleWeights(),
        NULL
    );


    particleCloud_.averagingM().setScalarSum
    (
        implicitEulerSource,
        partDatTmpImpl_,
        particleCloud_.particleWeights(),
        NULL
    );
    particleCloud_.makeSpecific(explicitEulerSource);
    particleCloud_.makeSpecific(implicitEulerSource);

    // limit source term
    scalar explicitEulerSourceInCell;
    forAll(explicitEulerSource,cellI)
    {
        explicitEulerSourceInCell = explicitEulerSource[cellI];

        if(mag(explicitEulerSourceInCell) > maxSource_ )
        {
             explicitEulerSource[cellI] = sign(explicitEulerSourceInCell) * maxSource_;
        }
    }

    //Reporting of integral quantities
    Field<scalar> writeValues; bool writeDiskNow=forceSubM(0).verboseToDisk(); //must call 'verboseToDisk()' only once since this function is incremeting a counter!
    writeValues.clear();
    if( forceSubM(0).verbose() || writeDiskNow)
    {
	    scalar dropletDepositionRate = gSum(-(explicitEulerSource
				                             +implicitEulerSource*fluidScalarField_)
			                                *explicitEulerSource.mesh().V()
					                       );
        writeValues.setSize(writeValues.size()+1, dropletDepositionRate);

        if(forceSubM(0).verbose())
          Info << "speciesID: " << speciesID 
               << ": total deposited amount of droplets [kg/s] (Eulerian) = " 
               <<  dropletDepositionRate << endl;

    }
    if( writeDiskNow )
      for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).verboseToDiskWrite(writeValues);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double scalarParticleFilter::LambdaParticleModel1(scalar Rep, scalar Stokes, double voidfraction) const
{
   if(!useDepositionRate_)  return 0;

   double partfract(0);
   double phip(0);
   double phip5(0);
   double funcPhi(0);
   double stokesCritical(0);
   double singlecollecEff(0);

   if(voidfraction>0.99) voidfraction = 0.99;
   if(voidfraction<0.10) voidfraction = 0.10;


   partfract      = 1.0 - voidfraction;
   phip           = pow(partfract, 0.333333333333333333);
   phip5          = phip*phip*phip*phip*phip;

   funcPhi        = (  6.0 - 6.0*phip5)
                  / (  6.0 - 9.0*phip 
                     + 9.0*phip5 
                     - 6.0*partfract*partfract 
                    );

   stokesCritical = Stokes * 0.5
                  * ( funcPhi+ 
                       1.14 
                     * pow(Rep,0.2) 
                     * pow(voidfraction,-1.5)
                    );

   double          StokesPowerThreePointTwo = pow(stokesCritical,3.2);

   singlecollecEff =  StokesPowerThreePointTwo 
                   /  (4.3 + StokesPowerThreePointTwo);

   return singlecollecEff;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double scalarParticleFilter::LambdaParticleModel2(scalar Rep, scalar Stokes, double voidfraction) const
{
   if(!useDepositionRate_)  return 0;

    //TODO: Implement an ALTERNATIVE filtration model on a per-particle basis
    double lambdaP(0);

    return lambdaP;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double scalarParticleFilter::stokesNumber(scalar magUr, scalar dParticle, scalar muFluid) const
{

    return   magUr
           * dropletDiameter_ * dropletDiameter_
           * dropletDensity_ 
           / (18.0 * dParticle * muFluid);
}

} // End namespace Foam

// ************************************************************************* //
