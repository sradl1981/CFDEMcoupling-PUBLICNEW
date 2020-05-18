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

    Two way general scalar exchange between DEM and CFD
    convective heat and species transfer model. The standard model is that of
    Deen, N.G. et al., Review of direct numerical simulation of 
    fluid–particle mass, momentum and heat transfer in dense gas–solid flows.
    Chemical Engineering Science 116 (2014) 710–724
    This correlation is based on that of Gunn (1978).

    Alternatively, the correclation of
    Li and Mason (2000), A computational investigation of transient heat
    transfer in  pneumatic transport of granular particles, Pow.Tech 112
    can be activated. However, this correlation is not suitable for 
    dense granular flows.

    This model allows an implicit/explicit split of the coupling term.
    The implicit/explicit splitting is realized in a force sub-model

    WARNING:
    This model REQUIRES the 'generalManual' speciesTransportModel

Class
    scalarGeneralExchangePhaseChange

SourceFiles
    scalarGeneralExchangePhaseChange.C

Contributing author and copyright holder of this model/file
    Copyright, 2015     Maryam Askarishahi, TU Graz (maryam.askarishahi@tugraz.at);
                        Stefan Radl, TU Graz (radl@tugraz.at)

\*---------------------------------------------------------------------------*/

#ifndef scalarGeneralExchangePhaseChange_H
#define scalarGeneralExchangePhaseChange_H

#include "scalarGeneralExchange.H"
#include "averagingModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scalarGeneralExchange Declaration
\*---------------------------------------------------------------------------*/

class scalarGeneralExchangePhaseChange
:
    public scalarGeneralExchange
{

private:

protected:

    //evaporation properties 
    scalarList          parameterVap_;      // parameters for vapor pressure
    scalar              Rvap_;  
    mutable bool        UseConstantEvaporationRate_;
    dimensionedScalar   deltaHEvap_;        //units: K; deltaHEvap = deltaH_v / cp_g   
    scalar              constantEvaporationRate_;      
    word                availableAreaModel_;    //Name of model for available area calculation
    scalar              LpStarNoEvaporation_;   //dimensionless liquid content below which there is no evaporation
    scalar              availableAreaLinearFactor_;  //dimensionless liquid content below which there is no evaporation
    scalar              dropletDensity_;
    scalar              dropletDiameter_;
    scalar              bParameter_;  // parameters for Stepanek Model 
    scalar              cParameter_;  

    //additional Lagrangian data arrays
    mutable double      **partTemp_;
    mutable double      **partDatFluxLatentHeat_;     // Lagrangian array CFDEM-->ExtCode for LATENT HEAT, 
                                                      // will be already registered by scalarGeneralExchange @ partHeatFluxPositionInRegister_

    void allocateMyArrays(scalar initialValue) const;

    void setPointersToExternalArrays(   word nameFlux,          int positionFlux,
                                        word nameTransCoeff,    int positionTransCoeff,
                                        word nameFluid,         int positionFluid
                                    ) const;

    inline double pVapor(double T) const
    {   //T in K, p in Pa, from http://ddbonline.ddbst.de/AntoineCalculation/AntoineCalculationCGI.exe?component=Water
            const double A = parameterVap_[0];
            const double B = parameterVap_[1];
            const double C = parameterVap_[2];
            const double Tmin = parameterVap_[3];
            const double Tmax = parameterVap_[4];
    
            double result = 0;
            if( T>Tmin )
                result = 133.322368 * pow(10.0, A - B / ( C + min(Tmax,T) ) );
            return result;
    }

    inline double expNegativeFast(double x) const
    {
        //approximates exp(-x), x MUST NOT include the minus!
        if(x > 2.5) //a large negative value must give 0
            return 0;
        else
            return 1.0-x*(1-x*0.5*(1-x*0.33333333*(1-x*0.25*(1-x*0.20*(1-x*0.166666667*(1-x*0.14286*(1-x*0.125)))))));
    }
    inline double logTaylorSeries(double x) const
    {
        return (x-1)*(1-(x-1)*(0.5+(x-1)*(0.333333-(x-1)*(0.25+(x-1)*(0.2-(x-1)*(0.166666667+(x-1)*(0.14286-(x-1)*0.125)))))));
    }

    mutable double (scalarGeneralExchangePhaseChange::*availableAreaFactor)(double LpStar, scalar dscaled) const;
    double availableAreaFactorLinear(double LpStar, scalar dscaled) const;
    double availabeAreaFactorStepanek(double LpStar, scalar dscaled) const;
    double availabeAreaFactorKariuki(double LpStar, scalar dscaled) const;


public:

    //- Runtime type information
    TypeName("scalarGeneralExchangePhaseChange");

    // Constructors

        //- Construct from components
        scalarGeneralExchangePhaseChange
        (
            const       dictionary& dict,
            cfdemCloud& sm
        );      

    // Destructor
        ~scalarGeneralExchangePhaseChange();


    // Member Functions
        void manipulateScalarField(volScalarField&, volScalarField&, int) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
