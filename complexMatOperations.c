#include "complexMatOperations.h"

#undef __FUNCT__
#define __FUNCT__ "complexProduct"
PetscErrorCode complexProduct(PetscScalar *resultReal, PetscScalar *resultImag, PetscScalar *factorsReal, PetscScalar *factorsImag, PetscInt nFactors){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscScalar partialResultReal, partialResultImag;
  PetscScalar partialFactorsReal[2], partialFactorsImag[2];

  /* The function is developed using recursivity */
  if(nFactors<1){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Number of factors must be a natural number"); //not valid number of factors
  }
  else if(nFactors==1){
    /* If only one factor, the result is such number */
    *resultReal = factorsReal[0];
    *resultImag = factorsImag[0];
    PetscFunctionReturn(0);
  }
  else if(nFactors==2){
    /* Base case of the recursion is when 2 factors */
    *resultReal = factorsReal[0]*factorsReal[1] - factorsImag[0]*factorsImag[1];
    *resultImag = factorsReal[0]*factorsImag[1] + factorsImag[0]*factorsReal[1];
    PetscFunctionReturn(0);
  }
  else{
    /* Product of elements from 1 to N-1 is calculated and stored in partialResultReal and partialResultImag */
    ierr = complexProduct(&partialResultReal, &partialResultImag, factorsReal, factorsImag, nFactors-1);CHKERRQ(ierr);

    /* Product is now calculated between partial results and the last element */
    partialFactorsReal[0]= partialResultReal;         partialFactorsImag[0]= partialResultImag;
    partialFactorsReal[1]= factorsReal[nFactors-1];   partialFactorsImag[1]= factorsImag[nFactors-1];
    ierr = complexProduct(resultReal, resultImag, partialFactorsReal, partialFactorsImag, 2);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
}

#undef __FUNCT__
#define __FUNCT__ "complexMatrixDeterminant_2x2"
PetscErrorCode complexMatrixDeterminant_2x2(PetscScalar *determinantReal, PetscScalar *determinantImag, PetscScalar *inputReal, PetscScalar *inputImag){
  PetscErrorCode ierr;
  PetscScalar auxReal1, auxReal2, auxImag1, auxImag2;
  PetscScalar factorsReal[2], factorsImag[2];

  PetscFunctionBegin;

  /* Product a11*a22 */
  factorsReal[0] = inputReal[0];  factorsImag[0] = inputImag[0];
  factorsReal[1] = inputReal[3];  factorsImag[1] = inputImag[3];
  ierr = complexProduct(&auxReal1,&auxImag1,factorsReal,factorsImag,2);CHKERRQ(ierr);

  /* Product a12*a21 */
  factorsReal[0] = inputReal[1];  factorsImag[0] = inputImag[1];
  factorsReal[1] = inputReal[2];  factorsImag[1] = inputImag[2];
  ierr = complexProduct(&auxReal2,&auxImag2,factorsReal,factorsImag,2);CHKERRQ(ierr);

  /* Det = a11*a22 - a12*a21 */
  *determinantReal = auxReal1 - auxReal2;
  *determinantImag = auxImag1 - auxImag2;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "complexMatrixDeterminant_3x3"
PetscErrorCode complexMatrixDeterminant_3x3(PetscScalar *determinantReal, PetscScalar *determinantImag, PetscScalar *inputReal, PetscScalar *inputImag){
  PetscErrorCode ierr;
  PetscInt i;
  PetscScalar factorsReal[2], factorsImag[2];
  PetscScalar subMatrixReal1[4], subMatrixReal2[4], subMatrixReal3[4], subMatrixImag1[4], subMatrixImag2[4], subMatrixImag3[4];
  PetscScalar subDetReal[3], subDetImag[3] ,auxReal[3], auxImag[3];

  PetscFunctionBegin;

  /* Sub-matrices assembly */

  subMatrixReal1[0] = inputReal[4]; subMatrixReal1[1] = inputReal[5];
  subMatrixReal1[2] = inputReal[7]; subMatrixReal1[3] = inputReal[8];
  subMatrixImag1[0] = inputImag[4]; subMatrixImag1[1] = inputImag[5];
  subMatrixImag1[2] = inputImag[7]; subMatrixImag1[3] = inputImag[8];

  subMatrixReal2[0] = inputReal[3]; subMatrixReal2[1] = inputReal[5];
  subMatrixReal2[2] = inputReal[6]; subMatrixReal2[3] = inputReal[8];
  subMatrixImag2[0] = inputImag[3]; subMatrixImag2[1] = inputImag[5];
  subMatrixImag2[2] = inputImag[6]; subMatrixImag2[3] = inputImag[8];

  subMatrixReal3[0] = inputReal[3]; subMatrixReal3[1] = inputReal[4];
  subMatrixReal3[2] = inputReal[6]; subMatrixReal3[3] = inputReal[7];
  subMatrixImag3[0] = inputImag[3]; subMatrixImag3[1] = inputImag[4];
  subMatrixImag3[2] = inputImag[6]; subMatrixImag3[3] = inputImag[7];

  /* Sub-matrices determinants */

  ierr = complexMatrixDeterminant_2x2(subDetReal, subDetImag, subMatrixReal1, subMatrixImag1);CHKERRQ(ierr);
  ierr = complexMatrixDeterminant_2x2(subDetReal+1, subDetImag+1, subMatrixReal2, subMatrixImag2);CHKERRQ(ierr);
  ierr = complexMatrixDeterminant_2x2(subDetReal+2, subDetImag+2, subMatrixReal3, subMatrixImag3);CHKERRQ(ierr);

  for(i=0; i<3; i++){
    /* Products a11*subdet1, a12*subdet2 and a13*subdet3 */
    factorsReal[0] = inputReal[i]; factorsImag[0] = inputImag[i];
    factorsReal[1] = subDetReal[i]; factorsImag[1] = subDetImag[i];
    ierr = complexProduct(auxReal+i,auxImag+i,factorsReal,factorsImag,2);CHKERRQ(ierr);
  }

  /* Determinant = a11*subdet1 - a12*subdet2 + a13*subdet3 */
  *determinantReal = auxReal[0] - auxReal[1] + auxReal[2];
  *determinantImag = auxImag[0] - auxImag[1] + auxImag[2];

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "complexScalarInverse"
PetscErrorCode complexScalarInverse(PetscScalar *resultReal, PetscScalar *resultImag, PetscScalar inputReal, PetscScalar inputImag){
  PetscFunctionBegin;

  if((inputReal==0) && (inputImag==0)){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Division by 0 is not possible");
  }
  *resultReal = inputReal/(inputReal*inputReal+inputImag*inputImag);
  *resultImag = -inputImag/(inputReal*inputReal+inputImag*inputImag);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "complexMatrixInverse_2x2"
PetscErrorCode complexMatrixInverse_2x2(PetscScalar *resultReal, PetscScalar *resultImag, PetscScalar *inputReal, PetscScalar *inputImag){
  PetscErrorCode ierr;
  PetscInt i;
  PetscScalar factorsReal[2], factorsImag[2];
  PetscScalar determinantReal, determinantImag, invDeterminantReal, invDeterminantImag;
  PetscScalar resultBufferReal[4], resultBufferImag[4];

  PetscFunctionBegin;

  /*Calculate the inverse of the determinant*/
  complexMatrixDeterminant_2x2(&determinantReal, &determinantImag, inputReal, inputImag);
  complexScalarInverse(&invDeterminantReal, &invDeterminantImag, determinantReal, determinantImag);

  /* Calulate product inv(det)*Mat */
  for(i=0; i<4; i++){
    /* Matrix elements (a11,a12, a21, a22 ) */
    factorsReal[0] = invDeterminantReal;  factorsImag[0] = invDeterminantImag;
    factorsReal[1] = inputReal[i];        factorsImag[1] = inputImag[i];
    ierr = complexProduct(resultBufferReal+i,resultBufferImag+i,factorsReal,factorsImag,2);CHKERRQ(ierr);
  }

  /*Reorder elements to get the proper expression of the inverse*/
  resultReal[0] = resultBufferReal[3];  resultReal[1] = -resultBufferReal[1];
  resultReal[2] = -resultBufferReal[2]; resultReal[3] = resultBufferReal[0];

  resultImag[0] = resultBufferImag[3];  resultImag[1] = -resultBufferImag[1];
  resultImag[2] = -resultBufferImag[2]; resultImag[3] = resultBufferImag[0];

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "complexMatrixInverse_3x3"
PetscErrorCode complexMatrixInverse_3x3(PetscScalar *resultReal, PetscScalar *resultImag, PetscScalar *inputReal, PetscScalar *inputImag){
  PetscErrorCode ierr;
  PetscInt i;
  PetscScalar factorsReal[2], factorsImag[2];
  PetscScalar auxMatrixReal[4], auxMatrixImag[4];
  PetscScalar determinantReal, determinantImag, invDeterminantReal, invDeterminantImag;
  PetscScalar resultBufferReal[9], resultBufferImag[9];

  PetscFunctionBegin;

  /* Calculate the inverse of the determinant */
  complexMatrixDeterminant_3x3(&determinantReal, &determinantImag, inputReal, inputImag);
  complexScalarInverse(&invDeterminantReal, &invDeterminantImag, determinantReal, determinantImag);

  /* Compute auxiliary matrix to multiply times the inverse of the determinant */

  /* a4*a8 - a5*a7 */
  auxMatrixReal[0] = inputReal[4];  auxMatrixReal[1] = inputReal[5];
  auxMatrixReal[2] = inputReal[7];  auxMatrixReal[3] = inputReal[8];
  auxMatrixImag[0] = inputImag[4];  auxMatrixImag[1] = inputImag[5];
  auxMatrixImag[2] = inputImag[7];  auxMatrixImag[3] = inputImag[8];
  ierr = complexMatrixDeterminant_2x2(resultBufferReal,resultBufferImag,auxMatrixReal,auxMatrixImag);CHKERRQ(ierr);

  /* a2*a7 - a1*a8 */
  auxMatrixReal[0] = inputReal[2];  auxMatrixReal[1] = inputReal[1];
  auxMatrixReal[2] = inputReal[8];  auxMatrixReal[3] = inputReal[7];
  auxMatrixImag[0] = inputImag[2];  auxMatrixImag[1] = inputImag[1];
  auxMatrixImag[2] = inputImag[8];  auxMatrixImag[3] = inputImag[7];
  ierr = complexMatrixDeterminant_2x2(resultBufferReal+1,resultBufferImag+1,auxMatrixReal,auxMatrixImag);CHKERRQ(ierr);

  /* a1*a5 - a2*a4 */
  auxMatrixReal[0] = inputReal[1];  auxMatrixReal[1] = inputReal[2];
  auxMatrixReal[2] = inputReal[4];  auxMatrixReal[3] = inputReal[5];
  auxMatrixImag[0] = inputImag[1];  auxMatrixImag[1] = inputImag[2];
  auxMatrixImag[2] = inputImag[4];  auxMatrixImag[3] = inputImag[5];
  ierr = complexMatrixDeterminant_2x2(resultBufferReal+2,resultBufferImag+2,auxMatrixReal,auxMatrixImag);CHKERRQ(ierr);

  /* a5*a6 - a3*a8 */
  auxMatrixReal[0] = inputReal[5];  auxMatrixReal[1] = inputReal[3];
  auxMatrixReal[2] = inputReal[8];  auxMatrixReal[3] = inputReal[6];
  auxMatrixImag[0] = inputImag[5];  auxMatrixImag[1] = inputImag[3];
  auxMatrixImag[2] = inputImag[8];  auxMatrixImag[3] = inputImag[6];
  ierr = complexMatrixDeterminant_2x2(resultBufferReal+3,resultBufferImag+3,auxMatrixReal,auxMatrixImag);CHKERRQ(ierr);

  /* a0*a8 - a2*a6 */
  auxMatrixReal[0] = inputReal[0];  auxMatrixReal[1] = inputReal[2];
  auxMatrixReal[2] = inputReal[6];  auxMatrixReal[3] = inputReal[8];
  auxMatrixImag[0] = inputImag[0];  auxMatrixImag[1] = inputImag[2];
  auxMatrixImag[2] = inputImag[6];  auxMatrixImag[3] = inputImag[8];
  ierr = complexMatrixDeterminant_2x2(resultBufferReal+4,resultBufferImag+4,auxMatrixReal,auxMatrixImag);CHKERRQ(ierr);

  /* a2*a3 - a0*a5 */
  auxMatrixReal[0] = inputReal[2];  auxMatrixReal[1] = inputReal[0];
  auxMatrixReal[2] = inputReal[5];  auxMatrixReal[3] = inputReal[3];
  auxMatrixImag[0] = inputImag[2];  auxMatrixImag[1] = inputImag[0];
  auxMatrixImag[2] = inputImag[5];  auxMatrixImag[3] = inputImag[3];
  ierr = complexMatrixDeterminant_2x2(resultBufferReal+5,resultBufferImag+5,auxMatrixReal,auxMatrixImag);CHKERRQ(ierr);

  /* a3*a7 - a4*a6 */
  auxMatrixReal[0] = inputReal[3];  auxMatrixReal[1] = inputReal[4];
  auxMatrixReal[2] = inputReal[6];  auxMatrixReal[3] = inputReal[7];
  auxMatrixImag[0] = inputImag[3];  auxMatrixImag[1] = inputImag[4];
  auxMatrixImag[2] = inputImag[6];  auxMatrixImag[3] = inputImag[7];
  ierr = complexMatrixDeterminant_2x2(resultBufferReal+6,resultBufferImag+6,auxMatrixReal,auxMatrixImag);CHKERRQ(ierr);

  /* a1*a6 - a0*a7 */
  auxMatrixReal[0] = inputReal[1];  auxMatrixReal[1] = inputReal[0];
  auxMatrixReal[2] = inputReal[7];  auxMatrixReal[3] = inputReal[6];
  auxMatrixImag[0] = inputImag[1];  auxMatrixImag[1] = inputImag[0];
  auxMatrixImag[2] = inputImag[7];  auxMatrixImag[3] = inputImag[6];
  ierr = complexMatrixDeterminant_2x2(resultBufferReal+7,resultBufferImag+7,auxMatrixReal,auxMatrixImag);CHKERRQ(ierr);

  /* a0*a4 - a1*a3 */
  auxMatrixReal[0] = inputReal[0];  auxMatrixReal[1] = inputReal[1];
  auxMatrixReal[2] = inputReal[3];  auxMatrixReal[3] = inputReal[4];
  auxMatrixImag[0] = inputImag[0];  auxMatrixImag[1] = inputImag[1];
  auxMatrixImag[2] = inputImag[3];  auxMatrixImag[3] = inputImag[4];
  ierr = complexMatrixDeterminant_2x2(resultBufferReal+8,resultBufferImag+8,auxMatrixReal,auxMatrixImag);CHKERRQ(ierr);

  /* (1/det) * AuxMatrix*/
  for(i=0; i<9; i++){
    factorsReal[0] = invDeterminantReal;  factorsImag[0] = invDeterminantImag;
    factorsReal[1] = resultBufferReal[i]; factorsImag[1] = resultBufferImag[i];
    ierr = complexProduct(resultReal+i,resultImag+i,factorsReal,factorsImag,2);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "zeroRowColIndexes_SymMat3x3"
PetscErrorCode zeroRowColIndexes_SymMat3x3(PetscInt *numZeros, PetscInt *indexes, PetscScalar *inputReal, PetscScalar *inputImag){
  indexes[0] = 0; indexes[1] = 0; indexes[2] = 0;
  *numZeros = 0;

  PetscFunctionBegin;

  /* Row/Col 0 is full of zeros */
  if((inputReal[0]==0)&&(inputReal[1]==0)&&(inputReal[2]==0)&&(inputImag[0]==0)&&(inputImag[1]==0)&&(inputImag[2]==0)){
    indexes[0] = 1; *numZeros+=1;
  }

  /* Row/Col 1 is full of zeros */
  if((inputReal[3]==0)&&(inputReal[4]==0)&&(inputReal[5]==0)&&(inputImag[3]==0)&&(inputImag[4]==0)&&(inputImag[5]==0)){
    indexes[1] = 1; *numZeros+=1;
  }

  /* Row/Col 2 is full of zeros */
  if((inputReal[6]==0)&&(inputReal[7]==0)&&(inputReal[8]==0)&&(inputImag[6]==0)&&(inputImag[7]==0)&&(inputImag[8]==0)){
    indexes[2] = 1; *numZeros+=1;
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "setUnion"
PetscErrorCode setUnion(PetscInt* output, PetscInt* input1, PetscInt sizeInput1, PetscInt* input2, PetscInt sizeInput2){
  PetscErrorCode ierr;
  PetscInt tmp,i=0,j=0,k=sizeInput1;
  PetscBool existent;

  PetscFunctionBegin;

  ierr = PetscMemcpy(output,input1,sizeInput1*sizeof(PetscInt));CHKERRQ(ierr);
  for (i=0; i<sizeInput2; i++){
    existent = PETSC_FALSE;
    for (j=0; j<sizeInput1; j++){
      if (input2[i] == input1[j])
      {
        existent =  PETSC_TRUE;
        break;
      }
    }
    if(!existent){
      output[k++]=input2[i];
    }
  }

  for (i=0; i<k; i++){
    for (j=i+1; j<k; j++){
      if (output[i] > output[j]){
        tmp =  output[i];
        output[i] = output[j];
        output[j] = tmp;
      }
    }
  }

  output[k]=-1;

  PetscFunctionReturn(0);
}
