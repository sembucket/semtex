/*
 * This file just includes one of the driver files
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include "pl/config.h"

#if defined(DRIVER_nil)
#  include "driver.nil"
#endif

#if defined(DRIVER_sm)
#  include "driver.sm"
#endif

#if defined(DRIVER_plot)
#  include "driver.plot"
#endif

