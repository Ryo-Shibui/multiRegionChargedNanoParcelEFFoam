/*---------------------------------------------------------------------------*\
  Maximum surface-film Courant number across all fluid regions.
\*---------------------------------------------------------------------------*/

scalar filmCoNum = -GREAT;
forAll(surfaceFilmFluid, regionI)
{
    filmCoNum = max(filmCoNum, surfaceFilmFluid[regionI].CourantNumber());
}

if (filmCoNum == -GREAT)
{
    filmCoNum = SMALL;
}
