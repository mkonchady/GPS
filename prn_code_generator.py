#*------------------------------------------------------------
#*- GPS PRN Code generator using Gold Codes
#*------------------------------------------------------------

#*-- Gold Codes for Satellite Vehicles
SV = {
   1: [2,6],
   2: [3,7],
   3: [4,8],
   4: [5,9],
   5: [1,9],
   6: [2,10],
   7: [1,8],
   8: [2,9],
   9: [3,10],
  10: [2,3],
  11: [3,4],
  12: [5,6],
  13: [6,7],
  14: [7,8],
  15: [8,9],
  16: [9,10],
  17: [1,4],
  18: [2,5],
  19: [3,6],
  20: [4,7],
  21: [5,8],
  22: [6,9],
  23: [1,3],
  24: [4,6],
  25: [5,7],
  26: [6,8],
  27: [7,9],
  28: [8,10],
  29: [1,6],
  30: [2,7],
  31: [3,8],
  32: [4,9],
}

#*---------------------------------------------------------------------------------------
#*- Parameters:
#*-    register: 10 element array
#*-    feedback: List of positions that are added modulo 2 and pre-pended to position 0
#*-    result:   List of positions that are added module 2 and returned to caller
#*-  Returns the bit that will be suffixed to the PRN code
#*----------------------------------------------------------------------------------------
def shift(register, feedback_pos, result_pos):    
    # calculate the value of the bit that will be added to the PRN code
    res_arr = [register[i-1] for i in result_pos]
    if len(res_arr) > 1:
        result = sum(res_arr) % 2
    else:
        result = res_arr[0]
        
    # calculate the feedback bit that will be inserted in position 0
    feedback = sum([register[i-1] for i in feedback_pos]) % 2
    
    # shift 1 position to the right
    for i in reversed(range(len(register[1:]))):
        register[i+1] = register[i]
    register[0] = feedback
    
    return result

#*---------------------------------------------------------------------------------------
#*- Parameters:
#*-    sv: Satellite Vehicle number
#*- Returns the PRN Code 
#*----------------------------------------------------------------------------------------
def PRN(sv):
    
    # initialize the two registers
    G1 = [1 for i in range(10)]
    G2 = [1 for i in range(10)]

    prn_code = [] # list of PRN code bits
    
    # Generate a bit that will be added to ca
    for i in range(1023):
        g1 = shift(G1, [3,10], [10])            # Generate the first part from the shift
        g2 = shift(G2, [2,3,6,8,9,10], SV[sv])  # Generate the second part from the shift for the particular satellite
        prn_code.append((g1 + g2) % 2)

    return prn_code


# Main 
print (PRN(24))