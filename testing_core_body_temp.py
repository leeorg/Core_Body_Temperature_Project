from core_body_temp import *

CBT_list = [36.6,36.64,36.67,36.7,36.75,36.79,36.82,36.83,36.84,36.88,36.95,37.03,37.07,37.1,37.12,
            37.14,37.16,37.18,37.2, 37.25,37.29,37.34,37.39,37.42,37.43,37.45,37.45,37.45,37.45,37.44,
            37.44,37.44,37.43,37.44,37.43,37.41,37.37,37.34,37.31,37.28,37.23,37.19,37.14,37.11,37.07,
            37.05,37.04,37.04,37.02,37.01]

def test_n_pt_mavg():
    assert n_pt_mavg(CBT_list, 21) == [36.769999999999996, 36.791666666666664, 36.81307692307692,
                                       36.833571428571425, 36.852666666666664, 36.870625000000004,
                                       36.887647058823532, 36.903888888888886, 36.91947368421053,
                                       36.936, 36.952857142857141, 36.988095238095241,
                                       37.023809523809526, 37.059523809523803, 37.094285714285711,
                                       37.127619047619042, 37.15904761904762, 37.189047619047621,
                                       37.218571428571437, 37.247142857142869, 37.273809523809533,
                                       37.297142857142866, 37.316190476190485, 37.333809523809521,
                                       37.349523809523802, 37.36333333333333, 37.374285714285712,
                                       37.382857142857141, 37.389047619047609, 37.392857142857146,
                                       37.391904761904762, 37.387142857142862, 37.377619047619042,
                                       37.364285714285721, 37.347619047619048, 37.329523809523806,
                                       37.310000000000002, 37.290476190476184, 37.269999999999996,
                                       37.249047619047609, 37.239499999999992, 37.228947368421046,
                                       37.217222222222219, 37.20470588235294, 37.189999999999998,
                                       37.173999999999999, 37.157142857142858, 37.14076923076923,
                                       37.124166666666667, 37.107272727272729]
    assert n_pt_mavg(CBT_list, 8) == [36.672000000000004, 36.69166666666667, 36.710000000000001,
                                     36.725000000000001, 36.737777777777779, 36.768888888888881,
                                     36.803333333333327, 36.843333333333334, 36.884444444444441,
                                     36.923333333333339, 36.960000000000001, 36.995555555555555,
                                     37.032222222222217, 37.07, 37.105555555555554,
                                     37.138888888888886, 37.167777777777786, 37.197777777777773,
                                     37.22999999999999, 37.263333333333335, 37.295555555555559,
                                     37.327777777777776, 37.357777777777777, 37.385555555555555,
                                     37.407777777777781, 37.424444444444447, 37.43555555555556,
                                     37.441111111111105, 37.44222222222222, 37.443333333333335,
                                     37.441111111111113, 37.43666666666666, 37.427777777777777,
                                     37.415555555555557, 37.401111111111113, 37.38333333333334,
                                     37.359999999999999, 37.333333333333336, 37.299999999999997,
                                     37.264444444444443, 37.226666666666667, 37.191111111111105,
                                     37.157777777777781, 37.12777777777778, 37.098888888888887,
                                     37.074444444444438, 37.060000000000002, 37.048571428571435,
                                     37.038333333333334, 37.031999999999996]
    print "passed n_pt_mavg"


def test_n_moving_stdev():
    assert n_moving_stdev(CBT_list, 21) == [0.10816653826391991, 0.12755272231393969, 0.14447517753260219,
                                            0.15858058581617804, 0.16976734780561381, 0.17905190122792106,
                                            0.18703373398148943, 0.19409333379130242, 0.20048479255739929,
                                            0.20866492908866774, 0.21755787407360902, 0.21747687342839553,
                                            0.21903141762030801, 0.21960137045283176, 0.21761368129259309,
                                            0.21584496199598377, 0.21224760921868516, 0.20637114047038582,
                                            0.19655242848374371, 0.18182801921439015, 0.16563442246333332,
                                            0.15166221866842342, 0.14119051634143856, 0.13169951368460703,
                                            0.1217159065396221, 0.11028750306962874, 0.097702171345954908,
                                            0.085037806718121126, 0.073478211866154886, 0.064741243202679821,
                                            0.067053638799027843, 0.077404318816385348, 0.094016209342047366,
                                            0.11056995200195374, 0.12692142301064729, 0.14090692639030214,
                                            0.15139352694220484, 0.15869707591749116, 0.1647118696390758,
                                            0.16860915639148272, 0.1670636500196695, 0.16465096878083577,
                                            0.16105554597370697, 0.15672897175772971, 0.14926486525636229,
                                            0.13958100566644849, 0.128028156243886, 0.11700865626536282,
                                            0.10500721475934507, 0.091442977761106117]
    assert n_moving_stdev(CBT_list, 8) == [0.057183913821982998, 0.070261416628663018,
                                              0.080415587212098114, 0.085690472882677976,
                                              0.08885068623507851, 0.083433273405225494,
                                              0.087464278422679773, 0.10024968827881742,
                                              0.10955718953029908, 0.11768602295939856,
                                              0.12227019260637426, 0.12299503151663324,
                                              0.11648795836670844, 0.1003742994994221,
                                              0.079074507761842258, 0.067720832179700041,
                                              0.070848037689440937, 0.080743076758596391,
                                              0.09367496997597706, 0.10259142264341681,
                                              0.10453601187044553, 0.10219806477837293,
                                              0.09257129384665902, 0.075184957124267274,
                                              0.057614620058146118, 0.037453675090406847,
                                              0.020069324297987738, 0.010540925533895495,
                                              0.0083333333333349916, 0.0070710678118674174,
                                              0.0078173595997071896, 0.01224744871391732,
                                              0.024381231397213803, 0.036438685181791199,
                                              0.049103066208852471, 0.060827625302980609,
                                              0.074999999999999525, 0.088459030064770572,
                                              0.099121138007994936, 0.103936412184459,
                                              0.10618380290797688, 0.10576441325470205,
                                              0.10009717500731272, 0.088568868370576259,
                                              0.073899330924650175, 0.060231036665308504,
                                              0.044721359549995975, 0.033380918415851051,
                                              0.021369760566432545, 0.016431676725154092]

    print "passed test_n_moving_stdev"
    
if __name__ == "__main__":
    print "*******************************************************"
    print "**** You are running testing_core_body_temp.py ********"
    print "*******************************************************"
    print
    test_n_pt_mavg()
    test_n_moving_stdev()
    print
    print
    print "YAY! All tests have passed!"
