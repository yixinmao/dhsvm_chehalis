/*
 * SUMMARY:      MassBalance.c - calculate basin-wide mass balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Mark Wigmosta
 * ORG:          Battelle - Pacific Northwest National Laboratory
 * E-MAIL:       ms_wigmosta@pnl.gov
 * ORIG-DATE:    Oct-96
 * DESCRIPTION:  Calculate water mass balance errors
 *               
 * DESCRIP-END.
 * FUNCTIONS:    MassBalance()
 * COMMENTS:
 * $Id: MassBalance.c,v 1.2 2002/10/01 21:30:59 nijssen Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "Calendar.h"

/*****************************************************************************
  Aggregate()
  
  Calculate the average values for the different fluxes and state variables
  over the basin.  Only the runoff is calculated as a total volume instead
  of an average.  In the current implementation the local radiation
  elements are not stored for the entire area.  Therefore these components
  are aggregated in AggregateRadiation() inside MassEnergyBalance().

  The aggregated values are set to zero in the function RestAggregate,
  which is executed at the beginning of each time step.
*****************************************************************************/
void MassBalance(DATE *Current, FILES *Out, AGGREGATED *Total,
		 WATERBALANCE *Mass)
{
  float NewWaterStorage;	/* water storage at the end of the time step */
  float Output;			/* total water flux leaving the basin;  */
  float MassError;		/* mass balance error m  */

  NewWaterStorage = Total->Runoff + Total->CanopyWater + Total->SoilWater +
    Total->Snow.Swq + Total->Soil.SatFlow;

  Output = Total->ChannelInt + Total->RoadInt + Total->Evap.ETot;

  MassError = (NewWaterStorage - Mass->OldWaterStorage) + Output -
    Total->Precip.Precip - Total->Snow.VaporMassFlux -
    Total->Snow.CanopyVaporMassFlux - Total->CulvertReturnFlow;

  /* update */
  Mass->OldWaterStorage = NewWaterStorage;
  Mass->CumPrecipIn += Total->Precip.Precip;
  Mass->CumRunoff += Total->Runoff;
  Mass->CumChannelInt += Total->ChannelInt;
  Mass->CumRoadInt += Total->RoadInt;
  Mass->CumET += Total->Evap.ETot;
  Mass->CumSnowVaporFlux += Total->Snow.VaporMassFlux +
    Total->Snow.CanopyVaporMassFlux;
  Mass->CumCulvertReturnFlow += Total->CulvertReturnFlow;
  Mass->CumCulvertToChannel += Total->CulvertToChannel;
  Mass->CumRunoffToChannel += Total->RunoffToChannel;

  PrintDate(Current, Out->FilePtr);

  fprintf(Out->FilePtr, " %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	  Total->Runoff, Total->CanopyWater, Total->SoilWater, Total->Snow.Swq,
	  Total->Soil.SatFlow, Total->ChannelInt, Total->RoadInt,
	  Total->CulvertReturnFlow, Total->Evap.ETot, Total->Precip.Precip,
	  Total->Snow.VaporMassFlux, Total->Snow.CanopyVaporMassFlux,
	  Mass->OldWaterStorage, Total->CulvertToChannel,
	  Total->RunoffToChannel, MassError);
}
