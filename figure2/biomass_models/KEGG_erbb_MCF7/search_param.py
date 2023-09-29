import numpy as np

from biomass.estimation import convert_scale, initialize_search_param

from .name2idx import C, V
from .ode import initial_values, param_values


class SearchParam(object):
    """Specify model parameters and/or initial values to optimize."""

    def __init__(self):
        # parameters
        self.idx_params = [
            # C.V1,
            # C.K1,
            C.V2,
            C.K2,
            # C.V3,
            # C.K3,
            C.V4,
            C.K4,
            # C.V5,
            # C.K5,
            C.V6,
            C.K6,
            # C.V7,
            # C.K7,
            C.V8,
            C.K8,
            # C.V9,
            # C.K9,
            # C.V10,
            # C.K10,
            C.V11,
            C.K11,
            # C.V12,
            # C.K12,
            # C.V13,
            # C.K13,
            # C.V14,
            # C.K14,
            # C.V15,
            # C.K15,
            C.V16,
            C.K16,
            C.V17,
            C.K17,
            C.V18,
            C.K18,
            C.V19,
            C.K19,
            # C.V20,
            # C.K20,
            # C.V21,
            # C.K21,
            # C.V22,
            # C.K22,
            # C.V23,
            # C.K23,
            # C.V24,
            # C.K24,
            # C.V25,
            # C.K25,
            # C.V26,
            # C.K26,
            # C.V27,
            # C.K27,
            # C.V28,
            # C.K28,
            C.V29,
            C.K29,
            C.V30,
            C.K30,
            C.V31,
            C.K31,
            C.V32,
            C.K32,
            C.V33,
            C.K33,
            # C.V34,
            # C.K34,
            # C.V35,
            # C.K35,
            # C.V36,
            # C.K36,
            # C.V37,
            # C.K37,
            # C.V38,
            # C.K38,
            # C.V39,
            # C.K39,
            # C.V40,
            # C.K40,
            # C.V41,
            # C.K41,
            # C.V42,
            # C.K42,
            # C.V43,
            # C.K43,
            # C.V44,
            # C.K44,
            # C.V45,
            # C.K45,
            C.V46,
            C.K46,
            C.V47,
            C.K47,
            C.V48,
            C.K48,
            C.V49,
            C.K49,
            C.V50,
            C.K50,
            C.V51,
            C.K51,
            # C.V52,
            # C.K52,
            C.V53,
            C.K53,
            C.V54,
            C.K54,
            C.V55,
            C.K55,
            # C.V56,
            # C.K56,
            C.V57,
            C.K57,
            C.V58,
            C.K58,
            C.V59,
            C.K59,
            C.V60,
            C.K60,
            C.V61,
            C.K61,
            C.V62,
            C.K62,
            C.V63,
            C.K63,
            C.V64,
            C.K64,
            C.V65,
            C.K65,
            C.V66,
            C.K66,
            C.V67,
            C.K67,
            # C.V68,
            # C.K68,
            # C.V69,
            # C.K69,
            # C.V70,
            # C.K70,
            # C.V71,
            # C.K71,
            C.V72,
            C.K72,
            # C.V73,
            # C.K73,
            # C.V74,
            # C.K74,
            # C.V75,
            # C.K75,
            # C.V76,
            # C.K76,
            # C.V77,
            # C.K77,
            # C.V78,
            # C.K78,
            # C.V79,
            # C.K79,
            # C.V80,
            # C.K80,
            # C.V81,
            # C.K81,
            # C.V82,
            # C.K82,
            # C.V83,
            # C.K83,
            # C.V84,
            # C.K84,
            C.V85,
            C.K85,
            C.V86,
            C.K86,
            # C.V87,
            # C.K87,
            # C.V88,
            # C.K88,
            # C.V89,
            # C.K89,
            # C.V90,
            # C.K90,
            C.V91,
            C.K91,
            C.V92,
            C.K92,
            # C.V93,
            # C.K93,
            # C.V94,
            # C.K94,
            # C.V95,
            # C.K95,
            # C.V96,
            # C.K96,
            # C.V97,
            # C.K97,
            # C.V98,
            # C.K98,
            C.V99,
            C.K99,
            C.V100,
            C.K100,
            # C.V101,
            # C.K101,
            # C.V102,
            # C.K102,
            C.V103,
            C.K103,
            C.V104,
            C.K104,
            # C.V105,
            # C.K105,
            # C.V106,
            # C.K106,
            # C.V107,
            # C.K107,
            # C.V108,
            # C.K108,
            C.V109,
            C.K109,
            C.V110,
            C.K110,
            # C.V111,
            # C.K111,
            C.kf112,
            C.kf113,
            C.kf114,
            C.kf115,
        ]

        # initial values
        self.idx_initials = [
            V.EGFR_EGFR,
            V.ERBB4_ERBB4,
            V.ERBB4_ERBB2,
            V.ERBB3_ERBB3,
            V.ERBB3_ERBB2,
            V.SHC2_SHC4_SHC3_SHC1,
            V.GRB2,
            V.PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3,
            V.SOS1_SOS2,
            V.GAB1,
            V.HRAS_KRAS_NRAS,
            V.ARAF_RAF1_BRAF,
            V.AKT3_AKT1_AKT2,
            V.MAP2K1_MAP2K2,
            V.MAPK1_MAPK3,
            V.MYC,
        ]

    def get_region(self):
        x = param_values()
        y0 = initial_values()

        search_rgn = np.zeros((2, len(x) + len(y0)))
        # Default: 0.1 ~ 10
        for i, j in enumerate(self.idx_params):
            if C.NAMES[j].startswith("K"):
                search_rgn[0, j] = 1e-2  # lower bound
                search_rgn[1, j] = 1e4 # upper bound
            elif C.NAMES[j].startswith("V"):
                search_rgn[0, j] = 1e-3  # lower bound
                search_rgn[1, j] = 1e3 # upper bound
            else:
                search_rgn[0, j] = 1e-5  # lower bound
                search_rgn[1, j] = 1 # upper bound
        # Default: 0.5 ~ 2
        for i, j in enumerate(self.idx_initials):
            search_rgn[0, j + len(x)] = 20  # lower bound
            search_rgn[1, j + len(x)] = 2000  # upper bound

        # search_rgn[:,C.parameter] = [lower_bound, upper_bound]
        # search_rgn[:,V.specie+len(x)] = [lower_bound, upper_bound]

        search_rgn = convert_scale(
            region=search_rgn,
            parameters=C.NAMES,
            species=V.NAMES,
            estimated_params=self.idx_params,
            estimated_initials=self.idx_initials,
        )

        return search_rgn

    def update(self, indiv):
        x = param_values()
        y0 = initial_values()

        for i, j in enumerate(self.idx_params):
            x[j] = indiv[i]
        for i, j in enumerate(self.idx_initials):
            y0[j] = indiv[i + len(self.idx_params)]

        # parameter constraints
        x[C.V3] = x[C.V1]
        x[C.K3] = x[C.K1]
        x[C.V5] = x[C.V1]
        x[C.K5] = x[C.K1]
        x[C.V7] = x[C.V1]
        x[C.K7] = x[C.K1]
        x[C.V9] = x[C.V1]
        x[C.K9] = x[C.K1]
        x[C.V10] = x[C.V1]
        x[C.K10] = x[C.K1]
        x[C.V12] = x[C.V1]
        x[C.K12] = x[C.K1]
        x[C.V13] = x[C.V1]
        x[C.K13] = x[C.K1]
        x[C.V14] = x[C.V1]
        x[C.K14] = x[C.K1]
        x[C.V15] = x[C.V1]
        x[C.K15] = x[C.K1]
        x[C.V20] = x[C.V1]
        x[C.K20] = x[C.K1]
        x[C.V21] = x[C.V1]
        x[C.K21] = x[C.K1]
        x[C.V22] = x[C.V1]
        x[C.K22] = x[C.K1]
        x[C.V23] = x[C.V1]
        x[C.K23] = x[C.K1]
        x[C.V24] = x[C.V1]
        x[C.K24] = x[C.K1]
        x[C.V25] = x[C.V1]
        x[C.K25] = x[C.K1]
        x[C.V26] = x[C.V1]
        x[C.K26] = x[C.K1]
        x[C.V27] = x[C.V1]
        x[C.K27] = x[C.K1]
        x[C.V28] = x[C.V1]
        x[C.K28] = x[C.K1]
        x[C.V34] = x[C.V1]
        x[C.K34] = x[C.K1]
        x[C.V35] = x[C.V2]
        x[C.K35] = x[C.K2]
        x[C.V36] = x[C.V1]
        x[C.K36] = x[C.K1]
        x[C.V37] = x[C.V2]
        x[C.K37] = x[C.K2]
        x[C.V38] = x[C.V1]
        x[C.K38] = x[C.K1]
        x[C.V39] = x[C.V2]
        x[C.K39] = x[C.K2]
        x[C.V40] = x[C.V1]
        x[C.K40] = x[C.K1]
        x[C.V41] = x[C.V2]
        x[C.K41] = x[C.K2]
        x[C.V42] = x[C.V1]
        x[C.K42] = x[C.K1]
        x[C.V43] = x[C.V2]
        x[C.K43] = x[C.K2]
        x[C.V44] = x[C.V1]
        x[C.K44] = x[C.K1]
        x[C.V45] = x[C.V2]
        x[C.K45] = x[C.K2]
        x[C.V52] = x[C.V1]
        x[C.K52] = x[C.K1]
        x[C.V56] = x[C.V1]
        x[C.K56] = x[C.K1]
        x[C.V68] = x[C.V1]
        x[C.K68] = x[C.K1]
        x[C.V69] = x[C.V2]
        x[C.K69] = x[C.K2]
        x[C.V70] = x[C.V1]
        x[C.K70] = x[C.K1]
        x[C.V71] = x[C.V2]
        x[C.K71] = x[C.K2]
        x[C.V73] = x[C.V1]
        x[C.K73] = x[C.K1]
        x[C.V74] = x[C.V2]
        x[C.K74] = x[C.K2]
        x[C.V75] = x[C.V1]
        x[C.K75] = x[C.K1]
        x[C.V76] = x[C.V2]
        x[C.K76] = x[C.K2]
        x[C.V77] = x[C.V1]
        x[C.K77] = x[C.K1]
        x[C.V78] = x[C.V2]
        x[C.K78] = x[C.K2]
        x[C.V79] = x[C.V1]
        x[C.K79] = x[C.K1]
        x[C.V80] = x[C.V2]
        x[C.K80] = x[C.K2]
        x[C.V81] = x[C.V1]
        x[C.K81] = x[C.K1]
        x[C.V82] = x[C.V2]
        x[C.K82] = x[C.K2]
        x[C.V83] = x[C.V1]
        x[C.K83] = x[C.K1]
        x[C.V84] = x[C.V2]
        x[C.K84] = x[C.K2]
        x[C.V87] = x[C.V1]
        x[C.K87] = x[C.K1]
        x[C.V88] = x[C.V2]
        x[C.K88] = x[C.K2]
        x[C.V89] = x[C.V1]
        x[C.K89] = x[C.K1]
        x[C.V90] = x[C.V2]
        x[C.K90] = x[C.K2]
        x[C.V93] = x[C.V1]
        x[C.K93] = x[C.K1]
        x[C.V94] = x[C.V2]
        x[C.K94] = x[C.K2]
        x[C.V95] = x[C.V1]
        x[C.K95] = x[C.K1]
        x[C.V96] = x[C.V2]
        x[C.K96] = x[C.K2]
        x[C.V97] = x[C.V1]
        x[C.K97] = x[C.K1]
        x[C.V98] = x[C.V2]
        x[C.K98] = x[C.K2]
        x[C.V101] = x[C.V1]
        x[C.K101] = x[C.K1]
        x[C.V102] = x[C.V2]
        x[C.K102] = x[C.K2]
        x[C.V105] = x[C.V1]
        x[C.K105] = x[C.K1]
        x[C.V106] = x[C.V2]
        x[C.K106] = x[C.K2]
        x[C.V107] = x[C.V1]
        x[C.K107] = x[C.K1]
        x[C.V108] = x[C.V2]
        x[C.K108] = x[C.K2]
        x[C.V111] = x[C.V1]
        x[C.K111] = x[C.K1]

        return x, y0
