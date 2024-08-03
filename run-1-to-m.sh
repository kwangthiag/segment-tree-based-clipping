# ------------------------------
# for ICPP
# 1-to-M
OMP_NUM_THREADS=10 bin/./seg data/ne_10m_ocean_0.txt data/continents_521.txt data/continents_1661.txt
OMP_NUM_THREADS=10 bin/./seg data/lakes_851348.txt data/parks_729699_-260.7_96.5_lakes_851348_fx.txt data/parks_844165_-230.6_13.8_lakes_851348_fx.txt data/parks_943452_-233.7_11.7_lakes_851348_fx.txt
OMP_NUM_THREADS=10 bin/./seg data/lakes_132372.txt  data/parks_169840_-232.4_13.1_lakes_132372_fx.txt data/parks_10613514_-235.8_9.1_lakes_132372_fx.txt data/parks_1243058_-225.3_15.4_lakes_132372_fx.txt
OMP_NUM_THREADS=10 bin/./seg data/lakes_174690.txt data/parks_321571_-226-8_-3-5_lakes_174690_fx.txt data/parks_34622_-230-6_-5-6_lakes_174690_fx.txt data/parks_140315_-234-2_-6-8_lakes_174690_fx.txt data/parks_169840_-230.2_-4.2_lakes_174690_fx.txt
# ------------------------------