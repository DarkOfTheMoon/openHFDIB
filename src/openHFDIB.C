/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________
                       | | | ||  ___|  _  \_   _| ___ \     H ybrid
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /     F ictitious
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \     D omain
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /     I mmersed
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/      B oundary
      | |
      |_|
-------------------------------------------------------------------------------
License

    openHFDIB is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Contributors
    Federico Municchi
\*---------------------------------------------------------------------------*/
#include "openHFDIB.H"
#include "fvMesh.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include <cmath>
#include <algorithm>

#include "interpolationCellPoint.H"
#include "interpolationCell.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
openHFDIB::openHFDIB(const Foam::fvMesh& mesh)
:
mesh_(mesh),
HFDIBDict_
 (
        IOobject
        (
            "HFDIBDict",
            "constant",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
 ),
 lambda_
 (
     IOobject
     (
         "lambda",
         mesh.time().timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("zero",dimless,0) // initialising to 0
 )
{
  //Initialize the library
  initialize();
}
//---------------------------------------------------------------------------//
openHFDIB::~openHFDIB()
{
 for(unsigned int i=0;i<immersedBodies_.size();i++)
  delete immersedBodies_[i];

 immersedBodies_.clear();
}
//---------------------------------------------------------------------------//
void openHFDIB::initialize()
{
 Info << "\nInitializing HFDIB library- reading dictionary\n";

 wordList stlNames( HFDIBDict_.lookup("stlNames") );

 HFDIBinterpDict_ = HFDIBDict_.subDict("interpolationSchemes");

 //Generate immersed objects
 Info << "\nInitializing HFDIB library- reading immersed objects\n";
 forAll(stlNames,name)
 {
  immersedBody * body_ptr = new immersedBody(stlNames[name],mesh_,HFDIBDict_);
  immersedBodies_.push_back(body_ptr);

 }

}
//---------------------------------------------------------------------------//
void openHFDIB::update()
{

Info << "\nHFDIB - updating immersed objects\n";

 for(unsigned int body_Id=0;body_Id<immersedBodies_.size();body_Id++)
  immersedBodies_[body_Id]->updateLambdaField(lambda_);

}
//---------------------------------------------------------------------------//
void openHFDIB::interpolateIB(volVectorField & V, volVectorField & Vs)
{

 //Create interpolator
 autoPtr<interpolation<vector> > interpV = interpolation<vector>::New(HFDIBinterpDict_, V);



 vector zeros = vector::zero;
 //Reset imposed field
 forAll(Vs,cellI)
 {
  Vs[cellI] =zeros;
 }

 //Loop over all the immersed bodies
 for(unsigned int body_Id=0;body_Id<immersedBodies_.size();body_Id++)
 {
  //Update imposed field according to lambda_
  immersedBodies_[body_Id]->updateVectoField(Vs, V.name());

  const std::vector<label> *                  surCells  = immersedBodies_[body_Id]->getSurfaceCellList();
  const std::vector< std::vector< point > > * intPoints = immersedBodies_[body_Id]->getInterpolationPoints();
  const std::vector< std::vector< label > > * intCells  = immersedBodies_[body_Id]->getInterpolationCells();

  //loop over all surface cells
  for(unsigned int scell=0;scell<surCells->size();scell++)
  {

   label cellI = (*surCells)[scell];
   //Check max order of accuracy
   bool allowedOrder[ORDER];

   for(int intPoint=0;intPoint<ORDER;intPoint++)
    if( (*intCells)[scell][intPoint] == -1 ) allowedOrder[intPoint] = false;
    else allowedOrder[intPoint] = true;

   bool  firstOrder = false;
   bool secondOrder = true;
   bool zeroOrder   = false;

   //Check is second order is possible
   if( allowedOrder[1] == false)
   {
    secondOrder = false;
    firstOrder  = true;
   }

   //Check if first order is possible
   if( allowedOrder[0] == false)
   {
    secondOrder = false;
    firstOrder  = false;
    zeroOrder   = true;
   }

   //Go for interpolation!
   if(secondOrder)
   {

     vector VP1 =  interpV->interpolate(  (*intPoints)[scell][1],
                                          (*intCells)[scell][0]
                                        ) - Vs[cellI];

     vector VP2 =  interpV->interpolate(  (*intPoints)[scell][2],
                                          (*intCells)[scell][1]
                                        ) - Vs[cellI];


    //distance between interpolation points
    double res_ = mag((*intPoints)[scell][2] -(*intPoints)[scell][1]);

    //cell center to surface distance
   double ds   = res_*(0.5-lambda_[cellI ]) ;

    vector quadCoeff = 1/(res_*res_) * ( VP2/2 - VP1 );
    vector linCoeff  = 1/(2*res_) * ( 4*VP1 - VP2 );
   // Info << "\nquadCoeff: " << quadCoeff << " linCoeff: " << linCoeff << " ds: " << ds << " res: " << res_;
    //Correct Vs
    Vs[cellI] = quadCoeff*ds*ds + linCoeff*ds + Vs[cellI]  ;
   }
   else if(firstOrder)
   {
     vector VP1 =  interpV->interpolate(  (*intPoints)[scell][1],
                                          (*intCells)[scell][0]
                                        ) - Vs[cellI];




    //distance between interpolation points
    double res_ = mag((*intPoints)[scell][1] -(*intPoints)[scell][0]);

    //cell center to surface distance
    double ds   = res_*(0.5-lambda_[cellI]) ;

    vector linCoeff = VP1/res_;
  //  Info << "\n linCoeff: " << linCoeff << " ds: " << ds << " res: " << res_;
   //Correct Vs
   Vs[cellI] = linCoeff*ds + Vs[cellI];
   }
   else if(zeroOrder)
   {
    //Zero order
    Vs[cellI] =Vs[cellI];
   }

  }

 }
}
//---------------------------------------------------------------------------//
void openHFDIB::updateIBForcing( volVectorField& V,
                      volVectorField& Vi,
                      fvVectorMatrix& VEqn,
                      volVectorField& IBforce,
                      volVectorField* explicitTerms
                    )
{

  //Reset force field
  IBforce = dimensionedVector("tmpIBsmallVec",IBforce.dimensions(),vector::zero);

  //Get diagonal coefficients of governing equations
  volScalarField rVA(VEqn.A());

  //Get RHS of governing equations
  volVectorField VH(VEqn.H());

  //Update boundary conditions
  rVA.correctBoundaryConditions();
  VH.correctBoundaryConditions();

  //Interpolate to evaluate the field value
  interpolateIB(V,Vi);

  //Loop over all the immersed bodies
  for(unsigned int body_Id=0;body_Id<immersedBodies_.size();body_Id++)
  {
   const std::vector<label> *   surCells  = immersedBodies_[body_Id]->getSurfaceCellList();

   //Force surface cells only
   for(unsigned int scell=0;scell<surCells->size();scell++)
   {
    label cellI = (*surCells)[scell];
    //Update IB forcing
     if(explicitTerms == NULL)
      IBforce[cellI] = rVA[cellI]*Vi[cellI] - VH[cellI];
     else
      IBforce[cellI] = rVA[cellI]*Vi[cellI] - VH[cellI] + (*explicitTerms)[cellI];
   }
  }


}
