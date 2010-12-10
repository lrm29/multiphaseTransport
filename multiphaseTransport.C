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

\*---------------------------------------------------------------------------*/

#include "multiphaseTransport.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//! Construct the multiphase::transport properties class.
//!
//! \param[in]  multiphaseTransportDictionary  Dictionary stored in constant/.
//! \param[in]  mesh    The mesh.
Foam::multiphase::transport::transport
(
    const IOdictionary& multiphaseTransportDictionary,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    multiphaseTransportDictionary_(multiphaseTransportDictionary),
    continuousPhaseSubDict_(multiphaseTransportDictionary.subDict("continuousPhase")),
    dispersedPhaseSubDict_(multiphaseTransportDictionary.subDict("dispersedPhase")),
    transportCoeffsSubDict_(multiphaseTransportDictionary.subDict("transportCoefficients")),
    dispersedPhases_(readLabel(multiphaseTransportDictionary_.lookup("dispersedPhases"))),
    Cl_(transportCoeffsSubDict_.lookup("Cl")),
    Cvm_(transportCoeffsSubDict_.lookup("Cvm")),
    rhoc_(continuousPhaseSubDict_.lookup("rho")),
    muc_(continuousPhaseSubDict_.lookup("mu")),
    nuc_(muc_/rhoc_),
    rhocField_(
            IOobject
            (
                    "rhoc",
                    mesh_.time().timeName(),
                    mesh_
             ),
             mesh_,
             dimensionedScalar(continuousPhaseSubDict_.lookup("rho"))
          ),
    mucField_(
            IOobject
            (
                    "muc",
                    mesh_.time().timeName(),
                    mesh_
             ),
             mesh_,
             dimensionedScalar(continuousPhaseSubDict_.lookup("mu"))
          ),
    nucField_(
            IOobject
            (
                    "nuc",
                    mesh_.time().timeName(),
                    mesh_
             ),
             mucField_/rhocField_
          )
{

    rhod_.setSize(dispersedPhases_);
    mud_.setSize(dispersedPhases_);
    nud_.setSize(dispersedPhases_);
    sigmad_.setSize(dispersedPhases_);
    rhodField_.setSize(dispersedPhases_);
    mudField_.setSize(dispersedPhases_);
    nudField_.setSize(dispersedPhases_);

    for (int i=0;i<dispersedPhases_;++i){
        rhod_.set(i, new dimensionedScalar(dispersedPhaseSubDict_.lookup("rho"+Foam::name(i+1))));
        mud_.set(i, new dimensionedScalar(dispersedPhaseSubDict_.lookup("mu"+Foam::name(i+1))));
        nud_.set(i, new dimensionedScalar(mud_[i]/rhod_[i]));
        sigmad_.set(i, new dimensionedScalar(dispersedPhaseSubDict_.lookup("sigma"+Foam::name(i+1))));
        rhodField_.set(i, new volScalarField(
                IOobject
                (
                    "rhodField"+ Foam::name(i+1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                 ),
                 mesh,
                 dimensionedScalar(dispersedPhaseSubDict_.lookup("rho"+Foam::name(i+1)))
             )
        );
        mudField_.set(i, new volScalarField(
                IOobject
                (
                    "mudField"+ Foam::name(i+1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                 ),
                 mesh,
                 dimensionedScalar(dispersedPhaseSubDict_.lookup("mu"+Foam::name(i+1)))
             )
        );
        nudField_.set(i, new volScalarField(
                IOobject
                (
                    "nudField"+ Foam::name(i+1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                 ),
                 mudField_[i]/rhodField_[i]
             )
        );
    }

}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
