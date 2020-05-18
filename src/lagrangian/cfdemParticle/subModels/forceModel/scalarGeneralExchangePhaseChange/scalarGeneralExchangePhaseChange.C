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

#include "scalarGeneralExchangePhaseChange.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"
#define ALARGECONCENTRATION 1e32

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(scalarGeneralExchangePhaseChange, 0);

addToRunTimeSelectionTable
(
    forceModel,
    scalarGeneralExchangePhaseChange,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
scalarGeneralExchangePhaseChange::scalarGeneralExchangePhaseChange
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    scalarGeneralExchange(dict,sm,"scalarGeneralExchangePhaseChange"),
    parameterVap_(propsDict_.lookup("parameterVap")),
    Rvap_(propsDict_.lookupOrDefault<scalar>("Rvap", 461.5)),
    UseConstantEvaporationRate_(false),
    deltaHEvap_(propsDict_.lookup("deltaHEvap")),
    constantEvaporationRate_(propsDict_.lookupOrDefault<scalar>("constantEvaporationRate", 0.0)),
    availableAreaModel_(propsDict_.lookupOrDefault<word>("availableAreaModel", "linear")),
    LpStarNoEvaporation_(propsDict_.lookupOrDefault<scalar>("LpStarNoEvaporation", 0.0)),
    availableAreaLinearFactor_(propsDict_.lookupOrDefault<scalar>("availableAreaLinearFactor", 1.0)),
    dropletDensity_(propsDict_.lookupOrDefault<scalar>("dropletDensity", 1000.0)),
    dropletDiameter_(propsDict_.lookupOrDefault<scalar>("dropletDiameter", 0)),
    bParameter_ (propsDict_.lookupOrDefault<scalar>("bParameter", 0.0)),
    cParameter_ (propsDict_.lookupOrDefault<scalar>("cParameter", 0,0)),
    partTemp_(NULL)
{
    //this model requires the partHeatFlux to be registered!
    if(partHeatFluxName_=="none")
        FatalError << "scalarGeneralExchangePhaseChange requires the 'partHeatFluxName' to be specified (i.e., not 'none') in order to push the latent heat flux. Set 'partHeatFluxName' in the dictionary correctly!" 
                   << abort(FatalError);

    if(availableAreaModel_=="linear")
        availableAreaFactor=&scalarGeneralExchangePhaseChange::availableAreaFactorLinear;
    else if(availableAreaModel_ == "stepanek")
        {
        availableAreaFactor = &scalarGeneralExchangePhaseChange::availabeAreaFactorStepanek;

        }
    else if(availableAreaModel_ == "kariuki")
    {
        availableAreaFactor = &scalarGeneralExchangePhaseChange::availabeAreaFactorKariuki;
        if(dropletDiameter_<=0)
            FatalError << "scalarGeneralExchangePhaseChange's availableAreaModel 'kariuki' requires dropletDiameter > 0! Please specify dropletDiameter > 0" 
                       << abort(FatalError);
        Info << "scalarGeneralExchangePhaseChange is using dropletDiameter: " << dropletDiameter_ << endl;

    }
    else
        FatalError << "scalarGeneralExchangePhaseChange: You have selected an invalid model for availableAreaModel. Check which models are available in the source code" 
                   << abort(FatalError);

    Info << "scalarGeneralExchangePhaseChange is using availableAreaModel: " << availableAreaModel_ << endl;
    Info << "scalarGeneralExchangePhaseChange is using LpStarNoEvaporation: " << LpStarNoEvaporation_ << endl;
    Info << "scalarGeneralExchangePhaseChange is using availableAreaLinearFactor: " << availableAreaLinearFactor_ << endl;

    allocateMyArrays(0.0);

    if (propsDict_.found("UseConstantEvaporationRate"))
    {
          UseConstantEvaporationRate_=readBool(propsDict_.lookup ("UseConstantEvaporationRate"));
          Info << "scalarGeneralExchangePhaseChange WARNING:  you are using a constant evaporation rate!" << endl;
    }

    particleCloud_.checkCG(true);

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

scalarGeneralExchangePhaseChange::~scalarGeneralExchangePhaseChange()
{
    delete partTemp_;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void scalarGeneralExchangePhaseChange::allocateMyArrays(scalar initVal) const
{
    scalarGeneralExchange::allocateMyArrays(initVal); //must call here!
    if(particleCloud_.numberOfParticlesChanged())
        particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1); 
}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void scalarGeneralExchangePhaseChange::manipulateScalarField(   volScalarField& explicitEulerSource,
                                                                volScalarField& implicitEulerSource,
                                                                int speciesID) const
{
    // reset Scalar field (== means hard reset)
    explicitEulerSource == dimensionedScalar("zero", explicitEulerSource.dimensions(), 0.);
    implicitEulerSource == dimensionedScalar("zero", implicitEulerSource.dimensions(), 0.);

    if(speciesID>=0 && particleSpeciesValue_[speciesID]<0.0)    //skip if species is not active
        return;

    if(speciesID<0)
        return; //no action in case here: latent heat source to be set during update of species fields

    //Set the names of the exchange fields
    word    fieldName;
    word    partDatName;
    scalar  transportParameter;

    fieldName          = eulerianFieldNames_[speciesID];
    partDatName        = partSpeciesNames_[speciesID];
    transportParameter = DMolecular_[speciesID];

    setPointersToExternalArrays( partSpeciesFluxNames_[speciesID],       partSpeciesFluxPositionInRegister_[speciesID],
                                 partSpeciesTransCoeffNames_[speciesID], partSpeciesTransCoeffPositionInRegister_[speciesID],
                                 partSpeciesFluidNames_[speciesID],      partSpeciesFluidPositionInRegister_[speciesID]
                               );
    if(probeIt_)
        particleCloud_.probeM().setOutputFile(typeName + "_" + fieldName + ".logDat");

    //==============================
    // get references
    const volScalarField& voidfraction_(particleCloud_.mesh().lookupObject<volScalarField> (voidfractionFieldName_));    // ref to voidfraction field
    const volVectorField& U_(particleCloud_.mesh().lookupObject<volVectorField> (velFieldName_));
    const volScalarField& fluidScalarField_(particleCloud_.mesh().lookupObject<volScalarField> (fieldName));            // ref to scalar field
    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField(); 
    //==============================

    // realloc the arrays and get data
    allocateMyArrays(0.0);

    // pull the temperature for each species
    particleCloud_.dataExchangeM().getData(partTempName_,"scalar-atom", partTemp_);

    if(particleSpeciesValue_[speciesID]>ALARGECONCENTRATION)   //only pull species if needed
          particleCloud_.dataExchangeM().getData(partDatName,"scalar-atom", partDat_);

    // calc La based heat flux
    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    scalar fluidValue(0);
    scalar fluidScalarFieldCell(0);
    label  cellI=0;
    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar dscaled(0);
    scalar dparcel(0);
    scalar numberParticlesInParcel(1);
    scalar nuf(0);
    scalar rhof(0);
    scalar magUr(0);
    scalar As(0);
    scalar Ad(0);
    scalar Lp(0);
    scalar LpStar(0);
    scalar availableArea(1); 
    scalar Rep(0);
    scalar Pr(0);

    double deltaT = voidfraction_.time().deltaTValue();

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    interpolationCellPoint<scalar> fluidScalarFieldInterpolator_(fluidScalarField_);

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                if(forceSubM(0).interpolation())
                {
                    position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = UInterpolator_.interpolate(position,cellI);
                    fluidValue = fluidScalarFieldInterpolator_.interpolate(position,cellI); //MUST use the mass loading, NOT partial density
                    fluidScalarFieldCell = fluidScalarField_[cellI];
                }else
                {
                    voidfraction = voidfraction_[cellI];
                    Ufluid       = U_[cellI];
                    fluidValue   = fluidScalarField_[cellI];
                    fluidScalarFieldCell = fluidValue;
                }

                // calc relative velocity
                Us      = particleCloud_.velocity(index);
                Ur      = Ufluid-Us;
                magUr   = mag(Ur);
                dscaled = 2*particleCloud_.radius(index);
                dparcel = dscaled;
                forceSubM(0).scaleDia(dscaled,index); //caution: this fct will scale ds!
                numberParticlesInParcel    = dparcel/dscaled;
                numberParticlesInParcel   *= numberParticlesInParcel*numberParticlesInParcel;
                As      = dscaled*dscaled*M_PI*numberParticlesInParcel;
                Lp      = max(0.0, partDat_[index][0]); //this is liquid VOLUME!
                LpStar  = Lp / (As*dscaled*0.166666666667); //divide by particle volume, 1/6 = 0.1666666666666667, since As already contains pi


                availableArea = (this->*availableAreaFactor)(LpStar, dscaled);
                Ad      = As * availableArea;
                nuf     = nufField[cellI];
                rhof    = rhoField[cellI];
                Rep     = dscaled*magUr/nuf;
                Pr      = max(SMALL,nuf/transportParameter); //This is Sc for species

                scalar Beta = transportParameter*(this->*Nusselt)(Rep,Pr,voidfraction)/(dscaled); 
             
                //MUST use the mass loading here, since this is how it is implemented in eulerianScalarField!
                scalar partMassLoadingSaturation = pVapor( partTemp_[index][0] ) /(partTemp_[index][0] * Rvap_) / rhof; 

                // calc transfer coefficient and rate
                scalar areaTimesTransferCoefficient = Beta 
                                                    * Ad;

                //Correct evaporation rate and transfer coefficients
                scalar tmpPartFlux = 0; //this is the MASS rate of evaporation!
                bool correctBeta = false;
                if(UseConstantEvaporationRate_)
                {
                    tmpPartFlux     =  -constantEvaporationRate_;
                    correctBeta     = true;
                }
		        else if ( tmpPartFlux < -Lp*dropletDensity_/deltaT )
		        {
                    //Bound particle flux on particles such that Lp cannot become negative
                    tmpPartFlux     = -Lp*dropletDensity_/deltaT;
                    correctBeta     = true;
		        }
		        else
		            tmpPartFlux   =  areaTimesTransferCoefficient * rhof //must consider fluid density here 
                                  * ( fluidValue - partMassLoadingSaturation );  //this is the MASS loading

                if(correctBeta)
		        {
                    if(mag(fluidValue - partMassLoadingSaturation)>0) //avoid division by zero
    			        areaTimesTransferCoefficient = tmpPartFlux / rhof
    			                                     / (fluidValue - partMassLoadingSaturation);
			        Beta = areaTimesTransferCoefficient
			             / max(Ad,1e-64); //avoid division by zero
		        }
                //split implicit/explicit contribution  
                //The source is in the transport equation for the MASS loading, so sources are DIVIDED by the fluid density
                forceSubM(0).explicitCorrScalar( partDatTmpImpl_[index][0], 
                                                 partDatTmpExpl_[index][0],
                                                 areaTimesTransferCoefficient,
                                                 fluidValue,
                                                 fluidScalarFieldCell,
                                                 partMassLoadingSaturation,  //the value on the particle DIVIDED by the is this!
                                                 forceSubM(0).verbose()
                                               );

                if(UseConstantEvaporationRate_)
                {
                    partDatTmpImpl_[index][0] = 0.0;
                    partDatTmpExpl_[index][0] = tmpPartFlux;
                }

                if(validPartFlux_)                                       
                    partDatFlux_[index][0]      += tmpPartFlux / dropletDensity_; //volumetric flux

                //Only for implicit handling on the DEM side!
                if(validPartTransCoeff_)
                    partDatTransCoeff_[index][0] = Beta / dropletDensity_; //volumetric flux
                if(validPartFluid_)
                    partDatFluid_[index][0]      = fluidValue; //WARNING: this is the mass loading!

                //add particle cooling due to evaporation from its surface
                partDatFluxLatentHeat_[index][0] = tmpPartFlux * deltaHEvap_.value(); //Pure explicit handling here!

                if( forceSubM(0).verbose())
                {
                    Pout << "fieldName = " << fieldName << endl;
                    Pout << "index    = " <<index << endl;
                    Pout << "LpStar   = " << LpStar << endl;
                    Pout << "availableArea = " << availableArea << endl;
                    Pout << "partFlux = " << tmpPartFlux << endl;
                    Pout << "magUr = " << magUr << endl;
                    Pout << "As = " << As << endl;
                    Pout << "r = " << particleCloud_.radius(index) << endl;
                    Pout << "dscaled = " << dscaled << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "Pr/Sc = " << Pr << endl;
                    Pout << "Nup/Shp = " << (this->*Nusselt)(Rep,Pr,voidfraction) << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "partDat_[index][0] = " << partDat_[index][0] << endl  ;
                    Pout << "partDatFluxLatentHeat_ = " << partDatFluxLatentHeat_[index][0] << endl  ;
                    Pout << "correctBeta            = " << correctBeta;
                    if(validPartFlux_)                                       
                        Pout << "partDatFlux_: " <<  partDatFlux_[index][0] << endl; 
                    if(validPartTransCoeff_)
                        Pout << "partDatTransCoeff: " <<  partDatTransCoeff_[index][0] << endl; 
                    Pout << "fluidValue = " << fluidValue << endl  ;
                }
                
                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.setSize(vValues.size()+1, Ur);
                    sValues.setSize(sValues.size()+1, (this->*Nusselt)(Rep,Pr,voidfraction));
                    sValues.setSize(sValues.size()+1, tmpPartFlux);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }
    }

    //Handle explicit and implicit source terms on the Euler side
    //these are simple summations!
    //MUST be mass-based!
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

    // scale with the cell volume to get (total) volume-specific source
    particleCloud_.makeSpecific(explicitEulerSource);
    explicitEulerSource*=-1;
    particleCloud_.makeSpecific(implicitEulerSource);
    implicitEulerSource*=-1;

    // limit explicit source term
    scalar explicitEulerSourceInCell;
    forAll(explicitEulerSource,cellI)
    {
        explicitEulerSourceInCell = explicitEulerSource[cellI]; //source per TOTAL volume!!
        if(mag(explicitEulerSourceInCell) > maxSource_ )
             explicitEulerSource[cellI] = sign(explicitEulerSourceInCell) * maxSource_;
    }


    //Reporting of integral quantities
    //these are NORMALIZED with the gas density!
    Field<scalar> writeValues; bool writeDiskNow=forceSubM(0).verboseToDisk(); //must call 'verboseToDisk()' only once since this function is incremeting a counter!
    writeValues.clear();
    if( forceSubM(0).verbose() || writeDiskNow)
    {
	    scalar exchangeRate = gSum(-(explicitEulerSource
                                        +implicitEulerSource*fluidScalarField_
                                     )
                                        *explicitEulerSource.mesh().V()
				                  );
        writeValues.setSize(writeValues.size()+1, exchangeRate);
        writeValues.setSize(writeValues.size()+1, exchangeRate * deltaHEvap_.value());

        if(forceSubM(0).verbose())
        {
            Info << "speciesID: " << speciesID 
                 << ": total phase-change particle-fluid species flux [kmol/s] (or [kg/s]) (Eulerian) = " 
                 << exchangeRate
                 << endl;
            Info << "total latent heat flux applied to particle cloud [W] (Eulerian) = " 
                 << exchangeRate * deltaHEvap_.value()
                 << endl;
        }
    }
    if( writeDiskNow )
      for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).verboseToDiskWrite(writeValues);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void scalarGeneralExchangePhaseChange::setPointersToExternalArrays( word nameFlux,          int positionFlux,
                                                                    word nameTransCoeff,    int positionTransCoeff,
                                                                    word nameFluid,         int positionFluid
                                                                  ) const
{
    //set standard pointers
    scalarGeneralExchange::setPointersToExternalArrays( nameFlux,       positionFlux,
                                                        nameTransCoeff, positionTransCoeff,
                                                        nameFluid,      positionFluid
                                                      );
    //set additional pointer for LATENT HEAT
    partDatFluxLatentHeat_  = particleCloud_.particleDatFieldsUserCFDEMToExt[partHeatFluxPositionInRegister_];

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double scalarGeneralExchangePhaseChange::availableAreaFactorLinear(double LpStar, scalar dscaled) const
{
    if(LpStar>LpStarNoEvaporation_)
        return min(1.,(LpStar-LpStarNoEvaporation_)*availableAreaLinearFactor_);
    else
        return 0.0;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double scalarGeneralExchangePhaseChange::availabeAreaFactorStepanek(double LpStar, scalar dscaled) const
{
    if(LpStar>LpStarNoEvaporation_)
    	return max(0., min(1., 1.0 / (1.0 + expNegativeFast( bParameter_ * (LpStar - cParameter_)))));
    else
        return 0.0;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double scalarGeneralExchangePhaseChange::availabeAreaFactorKariuki(double LpStar, scalar dscaled) const
{
    scalar VpToVd                   = (dscaled * dscaled *dscaled)/ (dropletDiameter_ * dropletDiameter_ * dropletDiameter_);
    scalar dropletFootPrintFraction = (dropletDiameter_ * dropletDiameter_ )/(4.0 * dscaled *dscaled);

    if(LpStar>LpStarNoEvaporation_)
        return max(0., min(1.,1.0 - expNegativeFast( -VpToVd * LpStar * logTaylorSeries(1-dropletFootPrintFraction))));
    else
        return 0.0;
}

} // End namespace Foam

// ************************************************************************* //
