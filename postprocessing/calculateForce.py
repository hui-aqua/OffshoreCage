BP1 = num_seg                       # Body attached point of fairlead 1 location
ML1a = BP1 - 1                      # 1st point on mooring line 1
ML1b = BP1 + num_seg                # 1st point on mooring line 2
BP2 = num_seg + (1*num_seg*2) + 1   # Body attached point of fairlead 2 location
ML2a = BP2 - 1                      # 1st point on mooring line 3
ML2b = BP2 + num_seg                # 1st point on mooring line 4
BP3 = num_seg + (2*num_seg*2) + 2   # Body attached point of fairlead 3 location
ML3a = BP3 - 1                      # 1st point on mooring line 5
ML3b = BP3 + num_seg                # 1st point on mooring line 6
BP4 = num_seg + (3*num_seg*2) + 3   # Body attached point of fairlead 4 location
ML4a = BP4 - 1                      # 1st point on mooring line 7
ML4b = BP4 + num_seg                # 1st point on mooring line 8

BP_14=[BP1,BP2,BP3,BP4]
ML_18=[ML1a,ML1b,ML2a,ML2b,ML3a,ML3b,ML4a,ML4b]  
