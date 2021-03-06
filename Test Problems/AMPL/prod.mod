# Model PROD from June 1989 version of CSTR 133


#     This model determines a series of workforce levels that will most
#economically meet demands and inventory requirements over time.  The
#formulation is motivated by the experiences of a large producer in the
#United States.  The data are for three products and 13 periods.


####  PRODUCTION SETS AND PARAMETERS  ###

set prd 'products';    # Members of the product group

param pt 'production time' {prd} > 0;
			# Crew-hours to produce 1000 units

param pc 'production cost' {prd} > 0;
			# Nominal production cost per 1000, used
			# to compute inventory and shortage costs


###  TIME PERIOD SETS AND PARAMETERS  ###

param first > 0 integer;
			# Index of first production period to be modeled

param last > first integer;
			# Index of last production period to be modeled

set time 'planning horizon' := first..last;


###  EMPLOYMENT PARAMETERS  ###

param cs 'crew size' > 0 integer;
			# Workers per crew

param sl 'shift length' > 0;
			# Regular-time hours per shift

param rtr 'regular time rate' > 0;
			# Wage per hour for regular-time labor

param otr 'overtime rate' > rtr;
			# Wage per hour for overtime labor

param iw 'initial workforce' >= 0 integer;
			# Crews employed at start of first period

param dpp 'days per period' {time} > 0;
			# Regular working days in a production period

param ol 'overtime limit' {time} >= 0;
			# Maximum crew-hours of overtime in a period

param cmin 'crew minimum' {time} >= 0;
			# Lower limit on average employment in a period

param cmax 'crew maximum' {t in time} >= cmin[t];
			# Upper limit on average employment in a period

param hc 'hiring cost' {time} >= 0;
			# Penalty cost of hiring a crew

param lc 'layoff cost' {time} >= 0;
			# Penalty cost of laying off a crew



###  DEMAND PARAMETERS  ###

param dem 'demand' {prd,first..last+1} >= 0;
			# Requirements (in 1000s)
			# to be met from current production and inventory

param pro 'promoted' {prd,first..last+1} logical;
			# true if product will be the subject
			# of a special promotion in the period


###  INVENTORY AND SHORTAGE PARAMETERS  ###

param rir 'regular inventory ratio' >= 0;
			# Proportion of non-promoted demand
			# that must be in inventory the previous period

param pir 'promotional inventory ratio' >= 0;
			# Proportion of promoted demand
			# that must be in inventory the previous period

param life 'inventory lifetime' > 0 integer;
			# Upper limit on number of periods that
			# any product may sit in inventory

param cri 'inventory cost ratio' {prd} > 0;
			# Inventory cost per 1000 units is
			# cri times nominal production cost

param crs 'shortage cost ratio' {prd} > 0;
			# Shortage cost per 1000 units is
			# crs times nominal production cost

param iinv 'initial inventory' {prd} >= 0;
			# Inventory at start of first period; age unknown

param iil 'initial inventory left' {p in prd, t in time}
	      := iinv[p] less sum {v in first..t} dem[p,v];

			# Initial inventory still available for allocation
			# at end of period t

param minv 'minimum inventory' {p in prd, t in time}
	      := dem[p,t+1] * (if pro[p,t+1] then pir else rir);

			# Lower limit on inventory at end of period t


###  VARIABLES  ###

var Crews{first-1..last} >= 0;
			# Average number of crews employed in each period

var Hire{time} >= 0;    # Crews hired from previous to current period

var Layoff{time} >= 0;  # Crews laid off from previous to current period

var Rprd 'regular production' {prd,time} >= 0;
			# Production using regular-time labor, in 1000s

var Oprd 'overtime production' {prd,time} >= 0;
			# Production using overtime labor, in 1000s

var Inv 'inventory' {prd,time,1..life} >= 0;
			# Inv[p,t,a] is the amount of product p that is
			# a periods old -- produced in period (t+1)-a --
			# and still in storage at the end of period t

var Short 'shortage' {prd,time} >= 0;
			# Accumulated unsatisfied demand at the end of period t



###  OBJECTIVE  ###

minimize cost:

    sum {t in time} rtr * sl * dpp[t] * cs * Crews[t] +
    sum {t in time} hc[t] * Hire[t] +
    sum {t in time} lc[t] * Layoff[t] +
    sum {t in time, p in prd} otr * cs * pt[p] * Oprd[p,t] +
    sum {t in time, p in prd, a in 1..life} cri[p] * pc[p] * Inv[p,t,a] +
    sum {t in time, p in prd} crs[p] * pc[p] * Short[p,t];

			# Full regular wages for all crews employed, plus
			# penalties for hiring and layoffs, plus
			# wages for any overtime worked, plus
			# inventory and shortage costs

			# (All other production costs are assumed
			# to depend on initial inventory and on demands,
			# and so are not included explicitly.)


###  CONSTRAINTS  ###

rlim 'regular-time limit' {t in time}:

    sum {p in prd} pt[p] * Rprd[p,t] <= sl * dpp[t] * Crews[t];

			# Hours needed to accomplish all regular-time
			# production in a period must not exceed
			# hours available on all shifts

olim 'overtime limit' {t in time}:

    sum {p in prd} pt[p] * Oprd[p,t] <= ol[t];

			# Hours needed to accomplish all overtime
			# production in a period must not exceed
			# the specified overtime limit

empl0 'initial crew level':  Crews[first-1] = iw;

			# Use given initial workforce

empl 'crew levels' {t in time}:  Crews[t] = Crews[t-1] + Hire[t] - Layoff[t];

			# Workforce changes by hiring or layoffs

emplbnd 'crew limits' {t in time}:  cmin[t] <= Crews[t] <= cmax[t];

			# Workforce must remain within specified bounds

dreq1 'first demand requirement' {p in prd}:

    Rprd[p,first] + Oprd[p,first] + Short[p,first]
			     - Inv[p,first,1] = dem[p,first] less iinv[p];

dreq 'demand requirements' {p in prd, t in first+1..last}:

    Rprd[p,t] + Oprd[p,t] + Short[p,t] - Short[p,t-1]
			  + sum {a in 1..life} (Inv[p,t-1,a] - Inv[p,t,a])
			                          = dem[p,t] less iil[p,t-1];

			# Production plus increase in shortage plus
			# decrease in inventory must equal demand

ireq 'inventory requirements' {p in prd, t in time}:

    sum {a in 1..life} Inv[p,t,a] + iil[p,t] >= minv[p,t];

			# Inventory in storage at end of period t
			# must meet specified minimum

izero 'impossible inventories' {p in prd, v in 1..life-1, a in v+1..life}:

    Inv[p,first+v-1,a] = 0;

			# In the vth period (starting from first)
			# no inventory may be more than v periods old
			# (initial inventories are handled separately)

ilim1 'new-inventory limits' {p in prd, t in time}:

    Inv[p,t,1] <= Rprd[p,t] + Oprd[p,t];

			# New inventory cannot exceed
			# production in the most recent period

ilim 'inventory limits' {p in prd, t in first+1..last, a in 2..life}:

    Inv[p,t,a] <= Inv[p,t-1,a-1];

			# Inventory left from period (t+1)-p
			# can only decrease as time goes on



