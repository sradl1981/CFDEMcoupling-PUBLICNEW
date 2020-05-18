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
    dense granular flows


Class
    scalarParticleFilter

SourceFiles
    scalarParticleFilter.C

Contributing author and copyright holder of this model/file
    Copyright, 2015     Stefan Radl, TU Graz (radl@tugraz.at)

Description
    This model is meant to implement interactions between CFDEM particles 
    and some Eulerian field variables (e.g., droplets suspended in the fluid)

    This model REQUIRES a suitable speciesTransportModel in order to draw
    the quantities of interest.
    
    The filtration coefficient is based on the work of XXXXXXXX (TODO!)

\*---------------------------------------------------------------------------*/

#ifndef scalarParticleFilter_H
#define scalarParticleFilter_H

#include "forceModel.H"
#include "averagingModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scalarParticleFilter Declaration
\*---------------------------------------------------------------------------*/

class scalarParticleFilter
:
    public forceModel
{
private:

    //Dictionaries to draw data from
    dictionary          propsDict_;

    const IOdictionary  scalarTransportProperties_;

    word                transportModelName_;
    dictionary          transportModelDict_;

    //Standard Eulerian field names - just to compute reference quantities
    word                voidfractionFieldName_;
    word                velFieldName_;

    mutable double      **partDat_;                   // Lagrangian array ExtCode-->CFDEM

    mutable double      **partDatFlux_;               // Lagrangian array CFDEM-->ExtCode

    mutable double      **partDatTmpExpl_;            // Lagrangian array - for explicit source to fluid eqns
    mutable double      **partDatTmpImpl_;            // Lagrangian array - for implicit source to fluid eqns
    
    
    //Switches
    bool                validPartFlux_;               //indicator if found, based on heat, but also used for species

    bool                useModel2_;                   //switch to activate calculation using Li-Mason
    bool                useDepositionRate_;           //switch to activate deposition calculation
    bool                useMomentumTransferRate_;     //switch to activate force calculation

    //Extra Fields - Will be manipulated / exchanged
    //Eulerian field names
    wordList            eulerianFieldNames_;
    //Lagrangian field names
    const List<bool>    isDropletCloud_;
    const wordList      partSpeciesNames_;
    const wordList      partSpeciesFluxNames_;

    //information related to external register
    int                 partHeatFluxPositionInRegister_;
    std::vector<int>    partSpeciesFluxPositionInRegister_;

    //Scalar properties
    scalar              maxSource_;                    // max (limited) value of src field

    //Droplet properties
    scalar              dropletDensity_;
    scalar              dropletDiameter_;

    void   allocateMyArrays() const;

      double LambdaParticleModel1(scalar Rep, scalar Stokes, double voidfraction) const;
      double LambdaParticleModel2(scalar Rep, scalar Stokes, double voidfraction) const;
      double (scalarParticleFilter::*LambdaParticle)(scalar Rep, scalar Stokes, double voidfraction) const;

      double stokesNumber(scalar magUr, scalar dParticle, scalar muFluid) const;

public:

    //- Runtime type information
    TypeName("scalarParticleFilter");

    // Constructors

        //- Construct from components
        scalarParticleFilter
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

          // Lagrangian array - for implicit source to fluid eqns
    
        ~scalarParticleFilter();


    // Member Functions
        void setForce() const;

        void manipulateScalarField(volScalarField&, volScalarField&, int) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
