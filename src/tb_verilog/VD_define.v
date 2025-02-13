
`define HAP_MAX_LENGTH              128
`define READ_MAX_LENGTH             64

`define HAP_LENGTH                  128
`define READ_LENGTH                 64

//* Score parameters
`define DP_PAIRHMM_SCORE_BITWIDTH   16
`define GENO_SRAM_WORD_AMOUNT       128

`define CONST_MATCH_SCORE           1
`define CONST_MISMATCH_SCORE        -4
`define CONST_GAP_OPEN              -6
`define CONST_GAP_EXTEND            -1


`define CONST_M2M                   -1
`define CONST_M2I                   -3072
`define CONST_I2I                   -1024
`define CONST_I2M                   -47

`define CONST_MATCH_BITWIDTH        16
`define CONST_BQ0_MATCH_SCORE       -443//-227
`define CONST_BQ1_MATCH_SCORE       -37//-18
`define CONST_BQ2_MATCH_SCORE       -1//-2
`define CONST_BQ3_MATCH_SCORE       0

`define CONST_MISMATCH_BITWIDTH     16
`define CONST_BQ0_MISMATCH_SCORE    -693//-896
`define CONST_BQ1_MISMATCH_SCORE    -1615//-1920
`define CONST_BQ2_MISMATCH_SCORE    -3043//-2944
`define CONST_BQ3_MISMATCH_SCORE    -4276//-3968