/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::multiphase::transport

Description
    Base-class for all transport models used by the incompressible turbulence
    models.

SourceFiles
    multiphaseTransport.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseTransport_H
#define multiphaseTransport_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//! \brief Namespace containing models to be used with multiphase solvers.
namespace multiphase
{

/*---------------------------------------------------------------------------*\
                Class multiphaseTransport Declaration
\*---------------------------------------------------------------------------*/

//! \defgroup multiphaseTransport Multiphase Transport Model.
//! This module contains documentation of the multiphase transport model.

//! \ingroup multiphaseTransport
//! \brief multiphaseTransport class
//!
//! multiphase::transport stores data on the physical properties of multiple
//! fluid phases so that they can be easily accessed and passed around.
//!
//! multiphaseTransportDictionary reads in values from the
//! multiphaseTransportProperties dictionary stored in the constant/ folder in the
//! case directory and should contain:
//! - dispersedPhases (number of dispersed phases)
//! - subDictionary continuousPhase
//!   - rho
//!   - mu
//! - subDictionary dispersedPhase (containing 'i' entries where i = dispersedPhases)
//!   - rhoi
//!   - mui
//!   - sigmai
//! - subDictionary transportCoefficients
//!   - Cvm (virtual mass coefficient)
//!   - Cl (lift coefficient)
//!
//! \todo I always use Newtonian transport models. If this ever changes then this
//! should be made into a base class or something. Look at the single phase
//! transport in OpenFOAM.
class transport
{

    // Private Data

        const fvMesh& mesh_;

        //! Dictionary containing changeable constants.
        const dictionary& multiphaseTransportDictionary_;
        const dictionary& continuousPhaseSubDict_;
        const dictionary& dispersedPhaseSubDict_;
        const dictionary& transportCoeffsSubDict_;

        const label dispersedPhases_;
        const dimensionedScalar Cl_;
        const dimensionedScalar Cvm_;

        dimensionedScalar rhoc_;
        PtrList<dimensionedScalar> rhod_;
        dimensionedScalar muc_;
        PtrList<dimensionedScalar> mud_;
        dimensionedScalar nuc_;
        PtrList<dimensionedScalar> nud_;

        PtrList<dimensionedScalar> sigmad_;

        volScalarField rhocField_;
        PtrList<volScalarField> rhodField_;
        volScalarField mucField_;
        PtrList<volScalarField> mudField_;
        volScalarField nucField_;
        PtrList<volScalarField> nudField_;


    // Private Member Functions

        //! Disallow copy construct
        transport(const transport&);

        //! Disallow default bitwise assignment
        void operator=(const transport&);


public:

    // Constructors

        //! Construct from components.
        transport
        (
            const IOdictionary& multiphaseTransportDictionary,
            const fvMesh& mesh
        );


    // Destructor
        ~transport(){};


    // Member Functions

        //! Return the density of the continuous phase.
        const dimensionedScalar& rhoc() const {return rhoc_;}
        //! Return a field for the density of the continuous phase.
        const volScalarField& rhocf() const {return rhocField_;}

        //! Return the density of the dispersed phase.
        const dimensionedScalar& rhod(const int& phaseNumber) const {return rhod_[phaseNumber];}
        //! Return a field for the density of the dispersed phase.
        const volScalarField& rhodf(const int& phaseNumber) const {return rhodField_[phaseNumber];}

        //! Return the dynamic viscosity of the continuous phase.
        const dimensionedScalar& muc() const {return muc_;}
        //! Return a field for the dynamic viscosity of the continuous phase.
        const volScalarField& mucf() const {return mucField_;}

        //! Return the dynamic viscosity of the dispersed phase.
        const dimensionedScalar& mud(const int& phaseNumber) const {return mud_[phaseNumber];}
        //! Return a field for the dynamic viscosity of the dispersed phase.
        const volScalarField& mudf(const int& phaseNumber) const {return mudField_[phaseNumber];}

        //! Return the kinematic viscosity of the continuous phase.
        const dimensionedScalar& nuc() const {return nuc_;}
        //! Return a field for the kinematic viscosity of the continuous phase.
        const volScalarField& nucf() const {return nucField_;}

        //! Return the kinematic viscosity of the dispersed phase.
        const dimensionedScalar& nud(const int& phaseNumber) const {return nud_[phaseNumber];}
        //! Return a field for the kinematic viscosity of the dispersed phase.
        const volScalarField& nudf(const int& phaseNumber) const {return nudField_[phaseNumber];}

        //! Return the interfacial tension of the ith dispersed phase with the continuous phase.
        const dimensionedScalar& sigmad(const int& phaseNumber) const {return sigmad_[phaseNumber];}

        //! Return the lift coefficient.
        const dimensionedScalar& Cl() const {return Cl_;}

        //! Return the virtual mass coefficient.
        const dimensionedScalar& Cvm() const {return Cvm_;}

        const label& Nd() const {return dispersedPhases_;}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace multiphase.

} // End namespace Foam.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
