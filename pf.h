#ifndef PF_H
#define PF_H

#include "complexMatOperations.h"
#include <petscsnes.h>
#include <petscdmnetwork.h>
#include <string.h>
#include <ctype.h>
#include <petscsys.h>
#include <petsctime.h>

#define MAXCHARSPERLINE 1000
#define PQ_BUS 1
#define PV_BUS 2
#define REF_BUS 3
#define ISOLATED_BUS 4

/* Bus data */
struct _p_VERTEXDATA{
  PetscScalar basekV; /* Bus Base kV */
  PetscInt    ide; /* Bus type code */
  PetscScalar gl; /* Active component of shunt admittance to ground */
  PetscScalar bl; /* Reactive component of shunt admittance to ground */
  PetscScalar vm; /* Bus voltage magnitude; in pu */
  PetscScalar va; /* Bus voltage phase angle */
  PetscInt    isConnBus;
  PetscInt    ngen; /* Number of generators incident at this bus */
  PetscInt    gidx[2]; /* list of inndices for accessing the generator data in GEN structure */
  PetscInt    nload;
  PetscInt    lidx[2];
  PetscInt    ndload;
  PetscInt    dlidx[2];
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _p_VERTEXDATA *VERTEXDATA;

/* Wye Load data */
struct _p_WYELOAD{
  PetscScalar 	pl; /* Active power component of constant MVA load */
  PetscScalar 	ql; /* Reactive power component of constant MVA load */
  PetscScalar 	ip; /* Active power component of constant current load: MW pu V */
  PetscScalar 	iq; /* Reactive power component of constant current load: Mvar pu V */
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _p_WYELOAD *WYELOAD;

/* Delta Load data */
struct _p_DELTALOAD{
  PetscScalar 	pl; /* Active power component of constant MVA load */
  PetscScalar 	ql; /* Reactive power component of constant MVA load */
  PetscScalar 	ip; /* Active power component of constant current load: MW pu V */
  PetscScalar 	iq; /* Reactive power component of constant current load: Mvar pu V */
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _p_DELTALOAD *DELTALOAD;

/* Generator data */
struct _p_GEN{
  PetscScalar 	pg; /* Generator active power output */
  PetscScalar 	qg; /* Generator reactive power output */
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _p_GEN *GEN;

/* Branch data */
struct _p_EDGEDATA{
  PetscInt 	status; /* Service status */
  PetscInt	internal_i; /* Internal From Bus Number */
  PetscInt	internal_j; /* Internal To Bus Number */
  PetscScalar   yff[2],yft[2],ytf[2],ytt[2]; /* [G,B] */
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _p_EDGEDATA *EDGEDATA;

/* PTI format data structure */
struct _p_PFDATA{
  PetscScalar sbase; /* System base MVA */
  PetscInt    nbus,ngen,nbranch,nwye,ndelta; /* # of buses,gens,branches, and loads */
  VERTEXDATA  bus;
  WYELOAD     wye;
  DELTALOAD   delta;
  GEN         gen;
  EDGEDATA    branch;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _p_PFDATA PFDATA;

PetscErrorCode DPFReadMatPowerData(PFDATA*, char*);

PetscErrorCode PrintDPF(PFDATA*);

#endif
