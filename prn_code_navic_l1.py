#*-------------------------------------------------------------------------------------------
#*- Script to generate PRN codes for NavIC in the L1 band
#*-  
#*- https://www.isro.gov.in/media_isro/pdf/SateliteNavigation/NavIC_SPS_ICD_L1_final.pdf
#*- 
#*- Page 30
#*- Test Data for L1 Band SPS IZ4 Data PRN codes for PRN IDs 1 to 16
#*-
#*- Columns: R0_Initial_Octal R1_Initial_Octal C_Initial_binary First_24_Chips Last_24_Chips
#*--------------------------------------------------------------------------------------------  
PRN_DATA = """
    0061727026503255544    0377627103341647600    10100    46555656    37436405
    1660130752435362260    0047555332635133703    10100    53142347    13443723
    0676457016477551225    0570574070736102152    00110    33721740    54403555
    1763467705267605701    0511013576745450615    10100    75664075    76447167
    1614265052776007236    1216243446624447775    10100    01536470    65327312
    1446113457553463523    0176452272675511054    00110    13424450    12250705
    1467417471470124574    0151055342317137706    10100    57666344    40632614
    0022513456555401603    1127720116046071664    00110    51073147    12652306
    0004420115402210365    0514407436155575524    00110    30526222    66712550
    0072276243316574510    0253070462740453542    00110    06645560    64762612
    1632356715721616750    0573371306324706336    10100    76773146    12063241
    1670164755420300763    1315135317732077306    00110    40543717    22556166
    1752127524253360255    1170303027726635012    10100    16630453    44473201
    0262220014044243135    1637171270537414673    00110    65441510    67166562
    1476157654546440020    0342370520251732111    00110    02211566    13320430
    1567545246612304745    0142423551056551362    10100    57305474    01701545
"""

lines = [line.strip() for line in PRN_DATA.strip().split('\n')] # split into lines and remove empty lines
PRN_DATA_TABLE = [line.split() for line in lines]               # create the table for testing
REG_SIZE = 55
REG_R0 = []
REG_R1 = []
REG_C = []

def octal_to_binary(octal_str):
    decimal_num = int(octal_str, 8)     # Convert octal string to decimal integer using base 8
    bin_num = bin(decimal_num)[2:]      # remove prefix 0b with 2:
    return list("0" * (REG_SIZE - len(bin_num)) + bin_num)    # pad with zeroes to reach len 55

def binary_to_octal(binary_str):
    octal_str = ""			# convert a binary string which must be a multiple of 3 to octal string
    if (len(binary_str) % 3 != 0):
        return -1
    for i in range(0, len(binary_str), 3):
        octal_value = int(binary_str[i]) * 4 + int(binary_str[i+1]) * 2 + int(binary_str[i+2]) * 1
        octal_str = octal_str + str(octal_value)
    return octal_str

def R0(index):
    global REG_R0
    return int(REG_R0[54- index])       # index from right to left

def set_R0(index, value):
    global REG_R0
    REG_R0[54 - index] = str(value)     # index from right to left

def R1(index):
    global REG_R1
    return int(REG_R1[54- index])       # index from right to left

def set_R1(index, value):
    global REG_R1
    REG_R1[54 - index] = str(value)     # index from right to left

def C(index):
    global REG_C
    return int(REG_C[4- index])         # index from right to left

def set_C(index, value):
    global REG_C
    REG_C[4 - index] = str(value)       # index from right to left

def gen_navic_sequence(prn_code):
    global REG_R0
    global REG_R1
    global REG_C
    global PRN_DATA_TABLE

    #*-- set the initial values of the registers and C
    R0_INIT = PRN_DATA_TABLE[prn_code-1][0]
    R1_INIT = PRN_DATA_TABLE[prn_code-1][1]
    C_INIT = PRN_DATA_TABLE[prn_code-1][2]

    #*-- initialize the R0 and R1 registers with the octal values and reverse order
    REG_R0 = octal_to_binary(R0_INIT)[::-1]
    REG_R1 = octal_to_binary(R1_INIT)[::-1]
    REG_C = list((C_INIT)[::-1])
    SEQUENCE = []
    PRN_CODE_LEN = 10230

    for _ in range(0,PRN_CODE_LEN):

        #*--- feedback for R0: Shift to the right by 1 and add feedback
        feedback = str(R0(50) ^ R0(45) ^ R0(40) ^ R0(20) ^ R0(10) ^ R0(5) ^ R0(0))
       
        #*--- calculate sigma 2 value
        sigma_2a = (R0(50) ^ R0(45) ^ R0(40)) & (R0(20) ^ R0(10) ^ R0(5) ^ R0(0))
        sigma_2b = ((R0(50) ^ R0(45)) & R0(40)) ^((R0(20) ^ R0(10)) & (R0(5) ^ R0(0)))
        sigma_2c = (R0(50) & R0(45)) ^ (R0(20) & R0(10)) ^ (R0(5) & R0(0))
        sigma2 = sigma_2a ^ sigma_2b ^ sigma_2c
    
        #*--- calculate r1a, shift register R0 and calculate r1b values
        r1a = sigma2 ^ (R0(40) ^ R0(35) ^ R0(30) ^ R0(25) ^ R0(15) ^ R0(0))
        REG_R0 = [feedback] + REG_R0[:-1]
        r1b = R1(50) ^ R1(45) ^ R1(40) ^ R1(20) ^ R1(10) ^ R1(5) ^ R1(0)

        #*--- feedback for R1: Shift to the right by 1 and add feedback
        feedback = r1a ^ r1b
        output = R1(0) ^ C(0)
        SEQUENCE.append(output)                 # append to the sequence
        REG_R1 = [str(feedback)] + REG_R1[:-1]  # shift R1 right and add feedback
        REG_C = [str(C(0))] + REG_C[:-1]        # circular shift C

    sequence_string = "".join(map(str, SEQUENCE))
    return binary_to_octal(sequence_string)

# ----------------------------------------
# Check all sequences in the table
# ----------------------------------------
if __name__ == "__main__":
    for i in range (1, 17):
        octal_str = gen_navic_sequence(i)
        first_24_chips = octal_str[:8]
        last_24_chips = octal_str[-8:]
        print ("Test: " + str(i))
        assert first_24_chips == PRN_DATA_TABLE[i - 1][3], "First 24 chips of " + str(i) + " did not match"
        assert last_24_chips == PRN_DATA_TABLE[i - 1][4], "Last 24 chips of " + str(i) + " did not match"
    print ("Unit tests passed")