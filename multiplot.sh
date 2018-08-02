

#!/bin/sh
pp-h5plot.py Ez_z_000000_0_3.h5 Ez_z_000000_p0_3.h5 \
Ez_z_000000_0_5.h5 Ez_z_000000_p0_5.h5 \
Ez_z_000000_0_6.h5 Ez_z_000000_p0_6.h5 \
Ez_z_000000_0_7.h5 Ez_z_000000_p0_7.h5 \
Ez_z_000000_1_0.h5 Ez_z_000000_p1_0.h5 \
--cno 0 0 1 1 5 5 7 7 2 2    \
--lstyle 0 1 0 1 0 1 0 1 0 1 \
--lab  "rb = 0.3" ""   "rb = 0.5" "" "rb = 0.6" "" "rb = 0.7" ""  "rb = 1.0" "" --latexoff 

#Ez_z_000000_1_2.h5  \
#Ez_z_000000_p1_2.h5 Ez_z_000000_1_5.h5 Ez_z_000000_p1_5.h5 \
#Ez_z_000000_2_0.h5 Ez_z_000000_p2_0.h5 \
#Ez_z_000000_0_3.h5 Ez_z_000000_p0_3.h5 Ez_z_000000_0_5.h5 \
#Ez_z_000000_p0_5.h5 Ez_z_000000_0_7.h5 Ez_z_000000_p0_7.h5 \
#Ez_z_000000_1_0.h5 Ez_z_000000_p1_0.h5 \
#Ez_z_000000_1_2.h5  \
#Ez_z_000000_p1_2.h5 Ez_z_000000_1_5.h5 Ez_z_000000_p1_5.h5 
#Ez_z_000000_2_0.h5 Ez_z_000000_p2_0.h5 \
#  "rb = 0.3" "" "rb = 0.5" "" "rb = 0.7" "" \
#"rb = 1.2" ""  "rb = 1.5" "" "rb = 2.0" ""
#"rb = 1.0" ""

# pp-h5plot.py low/ExmBy_x_000000_0_3.h5 low/ExmBy_x_000000_p0_3.h5 \
# low/ExmBy_x_000000_0_5.h5 low/ExmBy_x_000000_p0_5.h5 \
# low/ExmBy_x_000000_0_7.h5 low/ExmBy_x_000000_p0_7.h5 \
# low/ExmBy_x_000000_1_0.h5 low/ExmBy_x_000000_p1_0.h5 \
# --cno 0 0 1 1 5 5 7 7    \
# --lstyle 0 1 0 1 0 1 0 1 \
# --lab  "rb = 0.3" ""   "rb = 0.5" ""  "rb = 0.7" ""  "rb = 1.0" "" --latexoff 


# pp-h5plot.py high/ExmBy_x_000000_1_2.h5 high/ExmBy_x_000000_p1_2.h5 \
# high/ExmBy_x_000000_1_5.h5 high/ExmBy_x_000000_p1_5.h5 \
# high/ExmBy_x_000000_2_0.h5 high/ExmBy_x_000000_p2_0.h5 \
# high/ExmBy_x_000000_2_5.h5 high/ExmBy_x_000000_p2_5.h5 \