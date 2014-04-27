/*
 * SUMMARY:      ExecDump.c - Write selected output
 * USAGE:        Part of DHSVM
 *
 * DESCRIPTION:  Write selected output files
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIP-END.
 * FUNCTIONS:    ExecDump()
 *               DumpMaps()
 *               DumpPix()
 * COMMENTS:
 * $Id: ExecDump.c,v 1.2 2002/10/01 18:33:32 nijssen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "fileio.h"
#include "sizeofnt.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

/*****************************************************************************
  ExecDump()
*****************************************************************************/
void ExecDump(MAPSIZE * Map, DATE * Current, DATE * Start,
	      OPTIONSTRUCT * Options, DUMPSTRUCT * Dump, TOPOPIX ** TopoMap,
	      EVAPPIX ** EvapMap, PRECIPPIX ** PrecipMap,
	      RADCLASSPIX ** RadMap, SNOWPIX ** SnowMap, MET_MAP_PIX ** MetMap,
	      VEGPIX ** VegMap, LAYER * Veg, SOILPIX ** SoilMap, LAYER * Soil,
	      AGGREGATED * Total, UNITHYDRINFO * HydrographInfo,
	      Channel * ChannelData, float *Hydrograph)
{
  int i;			/* counter */
  int j;			/* counter */
  int x;
  int y;

  /* dump the aggregated basin values for this timestep */
  DumpPix(Current, IsEqualTime(Current, Start), &(Dump->Aggregate),
	  &(Total->Evap), &(Total->Precip), &(Total->RadClass), &(Total->Snow),
	  &(Total->Soil), Soil->MaxLayers, Veg->MaxLayers);
  fprintf(Dump->Aggregate.FilePtr, " %lu", Total->Saturated);
  fprintf(Dump->Aggregate.FilePtr, "\n");

  if (Options->Extent != POINT) {

    /* check whether the model state needs to be dumped at this timestep, and
       dump state if needed */

    if (Dump->NStates < 0) {
      StoreModelState(Dump->Path, Current, Map, Options, TopoMap, PrecipMap,
		      SnowMap, MetMap, RadMap, VegMap, Veg, SoilMap, Soil,
		      HydrographInfo, Hydrograph);
      if (Options->HasNetwork)
	StoreChannelState(Dump->Path, Current, ChannelData);
    }
    else {
      for (i = 0; i < Dump->NStates; i++) {
	if (IsEqualTime(Current, &(Dump->DState[i]))) {
	  StoreModelState(Dump->Path, Current, Map, Options, TopoMap,
			  PrecipMap, SnowMap, MetMap, RadMap, VegMap, Veg,
			  SoilMap, Soil, HydrographInfo, Hydrograph);
	  if (Options->HasNetwork)
	    StoreChannelState(Dump->Path, Current, ChannelData);
	}
      }
    }

    /* check which pixels need to be dumped, and dump if needed */

    for (i = 0; i < Dump->NPix; i++) {
      y = Dump->Pix[i].Loc.N;
      x = Dump->Pix[i].Loc.E;
      DumpPix(Current, IsEqualTime(Current, Start), &(Dump->Pix[i].OutFile),
	      &(EvapMap[y][x]), &(PrecipMap[y][x]), &(RadMap[y][x]),
	      &(SnowMap[y][x]), &(SoilMap[y][x]),
	      Soil->NLayers[(SoilMap[y][x].Soil - 1)],
	      Veg->NLayers[(VegMap[y][x].Veg - 1)]);

      fprintf(Dump->Pix[i].OutFile.FilePtr, "\n");
    }

    /* check which maps need to be dumped at this timestep, and dump maps 
       if needed */

    for (i = 0; i < Dump->NMaps; i++) {
      for (j = 0; j < Dump->DMap[i].N; j++) {
	if (IsEqualTime(Current, &(Dump->DMap[i].DumpDate[j]))) {
	  fprintf(stdout, "Dumping Maps at ");
	  PrintDate(Current, stdout);
	  fprintf(stdout, "\n");
	  DumpMap(Map, Current, &(Dump->DMap[i]), TopoMap, EvapMap,
		  PrecipMap, RadMap, SnowMap, SoilMap, Soil, VegMap, Veg);
	}
      }
    }
  }
}

/*****************************************************************************
  DumpMap()
*****************************************************************************/
void DumpMap(MAPSIZE * Map, DATE * Current, MAPDUMP * DMap, TOPOPIX ** TopoMap,
	     EVAPPIX ** EvapMap, PRECIPPIX ** PrecipMap, RADCLASSPIX ** RadMap,
	     SNOWPIX ** SnowMap, SOILPIX ** SoilMap, LAYER * Soil,
	     VEGPIX ** VegMap, LAYER * Veg)
{
  const char *Routine = "DumpMap";
  char DataLabel[MAXSTRING + 1];
  float Offset;
  float Range;
  int Index;
  int NSoil;			/* Number of soil layers for current pixel */
  int NVeg;			/* Number of veg layers for current pixel */
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  void *Array;

  sprintf(DataLabel, "%02d.%02d.%04d.%02d.%02d.%02d", Current->Month,
	  Current->Day, Current->Year, Current->Hour, Current->Min,
	  Current->Sec);

  /* find out what date we are dumping */
  for (Index = 0; Index < DMap->N; Index++) {
    if (IsEqualTime(Current, &(DMap->DumpDate[Index])))
      break;
  }

  switch (DMap->NumberType) {
  case NC_BYTE:
    if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_BYTE))))
      ReportError((char *) Routine, 1);
    break;
  case NC_CHAR:
    if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_CHAR))))
      ReportError((char *) Routine, 1);
    break;
  case NC_SHORT:
    if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_SHORT))))
      ReportError((char *) Routine, 1);
    break;
  case NC_INT:
    if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_INT))))
      ReportError((char *) Routine, 1);
    break;
    /* 8 bit integer not yet implemented in NetCDF 3.4, but anticipated in
       future versions */
    /*   case NC_LONG: */
    /*     if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_LONG))))
     */
    /*       ReportError((char *) Routine, 1); */
    /*     break; */
  case NC_FLOAT:
    if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_FLOAT))))
      ReportError((char *) Routine, 1);
    break;
  case NC_DOUBLE:
    if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_DOUBLE))))
      ReportError((char *) Routine, 1);
    break;
  default:
    Array = NULL;
    ReportError((char *) Routine, 40);
    break;
  }

  Offset = DMap->MinVal;
  Range = DMap->MaxVal - DMap->MinVal;

  switch (DMap->ID) {

  case 101:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].ETot;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((EvapMap[y][x].ETot - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 102:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      /* soil */
	      ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].EPot[NVeg];
	    else if (DMap->Layer <= NVeg)
	      /* vegetation layer */
	      ((float *) Array)[y * Map->NX + x] =
		EvapMap[y][x].EPot[DMap->Layer - 1];
	    else
	      /* vegetation layer not present at this pixel */
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EPot[NVeg] - Offset) /
				 Range * MAXUCHAR);
	    else if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EPot[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 103:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].EInt[NVeg];
	    else if (DMap->Layer <= NVeg)
	      ((float *) Array)[y * Map->NX + x] =
		EvapMap[y][x].EInt[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EInt[NVeg] - Offset) /
				 Range * MAXUCHAR);
	    else if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EInt[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 104:
    /* NETCDFWORK: This does not work for NETCDF.  Fix */
    if (DMap->Resolution == MAP_OUTPUT) {
      for (i = 0; i < Soil->MaxLayers; i++) {
	for (y = 0; y < Map->NY; y++) {
	  for (x = 0; x < Map->NX; x++) {
	    if (INBASIN(TopoMap[y][x].Mask)) {
	      NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	      if (DMap->Layer <= NVeg)
		((float *) Array)[y * Map->NX + x] =
		  EvapMap[y][x].ESoil[DMap->Layer - 1][i];
	      else
		((float *) Array)[y * Map->NX + x] = NA;
	    }
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	}
	Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY,
		      Map->NX, DMap, Index);
      }
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (i = 0; i < Soil->MaxLayers; i++) {
	for (y = 0; y < Map->NY; y++) {
	  for (x = 0; x < Map->NX; x++) {
	    if (INBASIN(TopoMap[y][x].Mask)) {
	      NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	      if (DMap->Layer <= NVeg)
		((unsigned char *) Array)[y * Map->NX + x] =
		  (unsigned char) ((EvapMap[y][x].ESoil[DMap->Layer - 1][i] -
				    Offset) / Range * MAXUCHAR);
	      else
		((unsigned char *) Array)[y * Map->NX + x] = 0;
	    }
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	}
	Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		      Index);
      }
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 105:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].EAct[NVeg];
	    else if (DMap->Layer <= NVeg)
	      ((float *) Array)[y * Map->NX + x] =
		EvapMap[y][x].EAct[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer > NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EAct[NVeg] - Offset) /
				 Range * MAXUCHAR);
	    else if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EAct[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 201:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = PrecipMap[y][x].Precip;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((PrecipMap[y][x].Precip - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 202:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer <= NVeg)
	      ((float *) Array)[y * Map->NX + x] =
		PrecipMap[y][x].IntRain[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((PrecipMap[y][x].IntRain[DMap->Layer - 1] -
				  Offset) / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 203:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer <= NVeg)
	      ((float *) Array)[y * Map->NX + x] =
		PrecipMap[y][x].IntSnow[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((PrecipMap[y][x].IntSnow[DMap->Layer - 1] -
				  Offset) / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 301:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = RadMap[y][x].Beam;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((RadMap[y][x].Beam - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 302:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = RadMap[y][x].Diffuse;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((RadMap[y][x].Diffuse - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 401:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] = SnowMap[y][x].HasSnow;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] = SnowMap[y][x].HasSnow;
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 402:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    SnowMap[y][x].SnowCoverOver;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    SnowMap[y][x].SnowCoverOver;
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 403:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned short *) Array)[y * Map->NX + x] = SnowMap[y][x].LastSnow;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) (((float) SnowMap[y][x].LastSnow - Offset) / Range
			     * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 404:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].Swq;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].Swq - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 405:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].Melt;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].Melt - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 406:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].PackWater;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].PackWater - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 407:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].TPack;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].TPack - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 408:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].SurfWater;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].SurfWater - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 409:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].TSurf;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].TSurf - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 410:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].ColdContent;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].ColdContent - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 501:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	    if (DMap->Layer <= NSoil)
	      ((float *) Array)[y * Map->NX + x] =
		SoilMap[y][x].Moist[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	    if (DMap->Layer <= NSoil)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((SoilMap[y][x].Moist[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 502:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	    if (DMap->Layer <= NSoil)
	      ((float *) Array)[y * Map->NX + x] =
		SoilMap[y][x].Perc[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	    if (DMap->Layer <= NSoil)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((SoilMap[y][x].Perc[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 503:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].TableDepth;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].TableDepth - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 504:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].SatFlow;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].SatFlow - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 505:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].TSurf;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].TSurf - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 506:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qnet;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qnet - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 507:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qs;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qs - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 508:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qe;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qe - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 509:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qg;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qg - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 510:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qst;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qst - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  default:
    ReportError((char *) Routine, 26);
    break;
  }

  free(Array);
}

/*****************************************************************************
  DumpPix()
*****************************************************************************/
void DumpPix(DATE * Current, int first, FILES * OutFile, EVAPPIX * Evap,
	     PRECIPPIX * Precip, RADCLASSPIX * Rad, SNOWPIX * Snow,
	     SOILPIX * Soil, int NSoil, int NVeg)
{
  int i;			/* counter */
  int j;			/* counter */

  if (first == 1) {
    fprintf(OutFile->FilePtr, "Date ");

    fprintf(OutFile->FilePtr,
	    "HasSnow OverSnow LastSnow Swq Melt PackWater TPack ");
    fprintf(OutFile->FilePtr, "SurfWater TSurf ColdContent ");
    fprintf(OutFile->FilePtr, "EvapTot ");
    for (i = 0; i < NVeg + 1; i++)
      fprintf(OutFile->FilePtr, "EPot%d ", i);
    for (i = 0; i < NVeg + 1; i++)
      fprintf(OutFile->FilePtr, "EAct%d ", i);
    for (i = 0; i < NVeg; i++)
      fprintf(OutFile->FilePtr, "EInt%d ", i);
    for (i = 0; i < NVeg; i++)
      for (j = 0; j < NSoil; j++)
	fprintf(OutFile->FilePtr, "ESoil%d%d ", i, j);
    fprintf(OutFile->FilePtr, "ESoil ");

    fprintf(OutFile->FilePtr, "Precip ");
    for (i = 0; i < NVeg; i++)
      fprintf(OutFile->FilePtr, "IntRain%d ", i);
    for (i = 0; i < NVeg; i++)
      fprintf(OutFile->FilePtr, "IntSnow%d ", i);

    fprintf(OutFile->FilePtr, "RadBeam RadDiff ");

    for (i = 0; i < NSoil; i++)
      fprintf(OutFile->FilePtr, "SoilMoist%d ", i);
    for (i = 0; i < NSoil; i++)
      fprintf(OutFile->FilePtr, "Perc%d ", i);
    fprintf(OutFile->FilePtr, "TableDepth SatFlow Runoff ");
    fprintf(OutFile->FilePtr, "SoilTemp Qnet Qs Qe Qg Qst Ra\n");
  }

  /* All variables are dumped in the case of a pixel dump */

  PrintDate(Current, OutFile->FilePtr);

  /* Snow */
  fprintf(OutFile->FilePtr, " %1d %1d %4d %g %g %g %g",
	  Snow->HasSnow, Snow->SnowCoverOver, Snow->LastSnow, Snow->Swq,
	  Snow->Melt, Snow->PackWater, Snow->TPack);
  fprintf(OutFile->FilePtr, " %g %g %g", Snow->SurfWater, Snow->TSurf,
	  Snow->ColdContent);
  fprintf(OutFile->FilePtr, " %g", Evap->ETot);
  for (i = 0; i < NVeg + 1; i++)
    fprintf(OutFile->FilePtr, " %g", Evap->EPot[i]);
  for (i = 0; i < NVeg + 1; i++)
    fprintf(OutFile->FilePtr, " %g", Evap->EAct[i]);
  for (i = 0; i < NVeg; i++)
    fprintf(OutFile->FilePtr, " %g", Evap->EInt[i]);
  for (i = 0; i < NVeg; i++)
    for (j = 0; j < NSoil; j++)
      fprintf(OutFile->FilePtr, " %g", Evap->ESoil[i][j]);
  fprintf(OutFile->FilePtr, " %g", Evap->EvapSoil);

  fprintf(OutFile->FilePtr, " %g", Precip->Precip);
  for (i = 0; i < NVeg; i++)
    fprintf(OutFile->FilePtr, " %g", Precip->IntRain[i]);
  for (i = 0; i < NVeg; i++)
    fprintf(OutFile->FilePtr, " %g", Precip->IntSnow[i]);

  fprintf(OutFile->FilePtr, " %g %g", Rad->Beam, Rad->Diffuse);

  for (i = 0; i < NSoil; i++)
    fprintf(OutFile->FilePtr, " %g", Soil->Moist[i]);
  for (i = 0; i < NSoil; i++)
    fprintf(OutFile->FilePtr, " %g", Soil->Perc[i]);
  fprintf(OutFile->FilePtr, " %g %g %g", Soil->TableDepth,
	  Soil->SatFlow, Soil->Runoff);
  fprintf(OutFile->FilePtr, " %g %g %g %g %g %g %g",
	  Soil->TSurf, Soil->Qnet, Soil->Qs, Soil->Qe, Soil->Qg, Soil->Qst,
	  Soil->Ra);
}
