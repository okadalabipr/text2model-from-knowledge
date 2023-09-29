from typing import Dict, List

from .name2idx import C, V


class ReactionNetwork(object):
    """
    Reaction indices grouped according to biological processes.
    This is used for sensitivity analysis (``target``='reaction').
    """

    def __init__(self) -> None:
        super(ReactionNetwork, self).__init__()
        self.reactions: Dict[str, List[int]] = {}

    @staticmethod
    def flux(t, y, x):
        """
        Flux vector.
        """

        v = {}
        v[1] = x[C.V1] * y[V.TGFA] * y[V.EGFR_EGFR] / (x[C.K1] + y[V.EGFR_EGFR])
        v[2] = x[C.V2] * y[V.a_EGFR_EGFR] / (x[C.K2] + y[V.a_EGFR_EGFR])
        v[3] = x[C.V3] * y[V.NRG3] * y[V.ERBB4_ERBB4] / (x[C.K3] + y[V.ERBB4_ERBB4])
        v[4] = x[C.V4] * y[V.a_ERBB4_ERBB4] / (x[C.K4] + y[V.a_ERBB4_ERBB4])
        v[5] = x[C.V5] * y[V.NRG3] * y[V.ERBB4_ERBB2] / (x[C.K5] + y[V.ERBB4_ERBB2])
        v[6] = x[C.V6] * y[V.a_ERBB4_ERBB2] / (x[C.K6] + y[V.a_ERBB4_ERBB2])
        v[7] = x[C.V7] * y[V.NRG2] * y[V.ERBB3_ERBB3] / (x[C.K7] + y[V.ERBB3_ERBB3])
        v[8] = x[C.V8] * y[V.a_ERBB3_ERBB3] / (x[C.K8] + y[V.a_ERBB3_ERBB3])
        v[9] = x[C.V9] * y[V.NRG2] * y[V.ERBB4_ERBB4] / (x[C.K9] + y[V.ERBB4_ERBB4])
        v[10] = x[C.V10] * y[V.NRG2] * y[V.ERBB3_ERBB2] / (x[C.K10] + y[V.ERBB3_ERBB2])
        v[11] = x[C.V11] * y[V.a_ERBB3_ERBB2] / (x[C.K11] + y[V.a_ERBB3_ERBB2])
        v[12] = x[C.V12] * y[V.NRG2] * y[V.ERBB4_ERBB2] / (x[C.K12] + y[V.ERBB4_ERBB2])
        v[13] = x[C.V13] * y[V.HBEGF] * y[V.EGFR_EGFR] / (x[C.K13] + y[V.EGFR_EGFR])
        v[14] = x[C.V14] * y[V.HBEGF] * y[V.ERBB4_ERBB4] / (x[C.K14] + y[V.ERBB4_ERBB4])
        v[15] = x[C.V15] * y[V.HBEGF] * y[V.ERBB4_ERBB2] / (x[C.K15] + y[V.ERBB4_ERBB2])
        v[16] = x[C.V16] * y[V.NRG1] * y[V.ERBB3_ERBB3] / (x[C.K16] + y[V.ERBB3_ERBB3])
        v[17] = x[C.V17] * y[V.NRG1] * y[V.ERBB4_ERBB4] / (x[C.K17] + y[V.ERBB4_ERBB4])
        v[18] = x[C.V18] * y[V.NRG1] * y[V.ERBB3_ERBB2] / (x[C.K18] + y[V.ERBB3_ERBB2])
        v[19] = x[C.V19] * y[V.NRG1] * y[V.ERBB4_ERBB2] / (x[C.K19] + y[V.ERBB4_ERBB2])
        v[20] = x[C.V20] * y[V.EREG] * y[V.EGFR_EGFR] / (x[C.K20] + y[V.EGFR_EGFR])
        v[21] = x[C.V21] * y[V.EREG] * y[V.ERBB4_ERBB4] / (x[C.K21] + y[V.ERBB4_ERBB4])
        v[22] = x[C.V22] * y[V.EREG] * y[V.ERBB4_ERBB2] / (x[C.K22] + y[V.ERBB4_ERBB2])
        v[23] = x[C.V23] * y[V.BTC] * y[V.EGFR_EGFR] / (x[C.K23] + y[V.EGFR_EGFR])
        v[24] = x[C.V24] * y[V.BTC] * y[V.ERBB4_ERBB4] / (x[C.K24] + y[V.ERBB4_ERBB4])
        v[25] = x[C.V25] * y[V.BTC] * y[V.ERBB4_ERBB2] / (x[C.K25] + y[V.ERBB4_ERBB2])
        v[26] = x[C.V26] * y[V.NRG4] * y[V.ERBB4_ERBB4] / (x[C.K26] + y[V.ERBB4_ERBB4])
        v[27] = x[C.V27] * y[V.NRG4] * y[V.ERBB4_ERBB2] / (x[C.K27] + y[V.ERBB4_ERBB2])
        v[28] = x[C.V28] * y[V.AREG] * y[V.EGFR_EGFR] / (x[C.K28] + y[V.EGFR_EGFR])
        v[29] = x[C.V29] * y[V.EGF] * y[V.EGFR_EGFR] / (x[C.K29] + y[V.EGFR_EGFR])
        v[30] = x[C.V30] * y[V.ERBB2_ERBB2] * y[V.SHC2_SHC4_SHC3_SHC1] / (x[C.K30] + y[V.SHC2_SHC4_SHC3_SHC1])
        v[31] = x[C.V31] * y[V.a_SHC2_SHC4_SHC3_SHC1] / (x[C.K31] + y[V.a_SHC2_SHC4_SHC3_SHC1])
        v[32] = x[C.V32] * y[V.ERBB2_ERBB2] * y[V.GRB2] / (x[C.K32] + y[V.GRB2])
        v[33] = x[C.V33] * y[V.a_GRB2] / (x[C.K33] + y[V.a_GRB2])
        v[34] = x[C.V34] * y[V.a_EGFR_EGFR] * y[V.PLCG1_PLCG2] / (x[C.K34] + y[V.PLCG1_PLCG2])
        v[35] = x[C.V35] * y[V.a_PLCG1_PLCG2] / (x[C.K35] + y[V.a_PLCG1_PLCG2])
        v[36] = x[C.V36] * y[V.a_EGFR_EGFR] * y[V.CBL_CBLB] / (x[C.K36] + y[V.CBL_CBLB])
        v[37] = x[C.V37] * y[V.a_CBL_CBLB] / (x[C.K37] + y[V.a_CBL_CBLB])
        v[38] = x[C.V38] * y[V.a_EGFR_EGFR] * y[V.STAT5A_STAT5B] / (x[C.K38] + y[V.STAT5A_STAT5B])
        v[39] = x[C.V39] * y[V.a_STAT5A_STAT5B] / (x[C.K39] + y[V.a_STAT5A_STAT5B])
        v[40] = x[C.V40] * y[V.a_EGFR_EGFR] * y[V.CRK_CRKL] / (x[C.K40] + y[V.CRK_CRKL])
        v[41] = x[C.V41] * y[V.a_CRK_CRKL] / (x[C.K41] + y[V.a_CRK_CRKL])
        v[42] = x[C.V42] * y[V.a_EGFR_EGFR] * y[V.SRC] / (x[C.K42] + y[V.SRC])
        v[43] = x[C.V43] * y[V.a_SRC] / (x[C.K43] + y[V.a_SRC])
        v[44] = x[C.V44] * y[V.a_EGFR_EGFR] * y[V.NCK1_NCK2] / (x[C.K44] + y[V.NCK1_NCK2])
        v[45] = x[C.V45] * y[V.a_NCK1_NCK2] / (x[C.K45] + y[V.a_NCK1_NCK2])
        v[46] = x[C.V46] * y[V.a_EGFR_EGFR] * y[V.GRB2] / (x[C.K46] + y[V.GRB2])
        v[47] = x[C.V47] * y[V.a_EGFR_EGFR] * y[V.SHC2_SHC4_SHC3_SHC1] / (x[C.K47] + y[V.SHC2_SHC4_SHC3_SHC1])
        v[48] = x[C.V48] * y[V.a_ERBB4_ERBB4] * y[V.SHC2_SHC4_SHC3_SHC1] / (x[C.K48] + y[V.SHC2_SHC4_SHC3_SHC1])
        v[49] = x[C.V49] * y[V.a_ERBB4_ERBB4] * y[V.GRB2] / (x[C.K49] + y[V.GRB2])
        v[50] = x[C.V50] * y[V.a_ERBB4_ERBB4] * y[V.PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3] / (x[C.K50] + y[V.PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3])
        v[51] = x[C.V51] * y[V.a_PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3] / (x[C.K51] + y[V.a_PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3])
        v[52] = x[C.V52] * y[V.a_ERBB4_ERBB4] * y[V.STAT5A_STAT5B] / (x[C.K52] + y[V.STAT5A_STAT5B])
        v[53] = x[C.V53] * y[V.a_ERBB4_ERBB2] * y[V.SHC2_SHC4_SHC3_SHC1] / (x[C.K53] + y[V.SHC2_SHC4_SHC3_SHC1])
        v[54] = x[C.V54] * y[V.a_ERBB4_ERBB2] * y[V.GRB2] / (x[C.K54] + y[V.GRB2])
        v[55] = x[C.V55] * y[V.a_ERBB4_ERBB2] * y[V.PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3] / (x[C.K55] + y[V.PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3])
        v[56] = x[C.V56] * y[V.a_ERBB4_ERBB2] * y[V.STAT5A_STAT5B] / (x[C.K56] + y[V.STAT5A_STAT5B])
        v[57] = x[C.V57] * y[V.a_ERBB3_ERBB3] * y[V.PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3] / (x[C.K57] + y[V.PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3])
        v[58] = x[C.V58] * y[V.a_ERBB3_ERBB2] * y[V.SHC2_SHC4_SHC3_SHC1] / (x[C.K58] + y[V.SHC2_SHC4_SHC3_SHC1])
        v[59] = x[C.V59] * y[V.a_ERBB3_ERBB2] * y[V.GRB2] / (x[C.K59] + y[V.GRB2])
        v[60] = x[C.V60] * y[V.a_ERBB3_ERBB2] * y[V.PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3] / (x[C.K60] + y[V.PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3])
        v[61] = x[C.V61] * y[V.a_SHC2_SHC4_SHC3_SHC1] * y[V.GRB2] / (x[C.K61] + y[V.GRB2])
        v[62] = x[C.V62] * y[V.a_PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3] * y[V.AKT3_AKT1_AKT2] / (x[C.K62] + y[V.AKT3_AKT1_AKT2])
        v[63] = x[C.V63] * y[V.a_AKT3_AKT1_AKT2] / (x[C.K63] + y[V.a_AKT3_AKT1_AKT2])
        v[64] = x[C.V64] * y[V.a_GRB2] * y[V.SOS1_SOS2] / (x[C.K64] + y[V.SOS1_SOS2])
        v[65] = x[C.V65] * y[V.a_SOS1_SOS2] / (x[C.K65] + y[V.a_SOS1_SOS2])
        v[66] = x[C.V66] * y[V.a_GRB2] * y[V.GAB1] / (x[C.K66] + y[V.GAB1])
        v[67] = x[C.V67] * y[V.a_GAB1] / (x[C.K67] + y[V.a_GAB1])
        v[68] = x[C.V68] * y[V.a_PLCG1_PLCG2] * y[V.CAMK2A_CAMK2B_CAMK2D_CAMK2G] / (x[C.K68] + y[V.CAMK2A_CAMK2B_CAMK2D_CAMK2G])
        v[69] = x[C.V69] * y[V.a_CAMK2A_CAMK2B_CAMK2D_CAMK2G] / (x[C.K69] + y[V.a_CAMK2A_CAMK2B_CAMK2D_CAMK2G])
        v[70] = x[C.V70] * y[V.a_PLCG1_PLCG2] * y[V.PRKCA_PRKCB_PRKCG] / (x[C.K70] + y[V.PRKCA_PRKCB_PRKCG])
        v[71] = x[C.V71] * y[V.a_PRKCA_PRKCB_PRKCG] / (x[C.K71] + y[V.a_PRKCA_PRKCB_PRKCG])
        v[72] = x[C.V72] * y[V.a_GAB1] * y[V.PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3] / (x[C.K72] + y[V.PIK3CA_PIK3CB_PIK3CD_PIK3R1_PIK3R2_PIK3R3])
        v[73] = x[C.V73] * y[V.a_AKT3_AKT1_AKT2] * y[V.MTOR] / (x[C.K73] + y[V.MTOR])
        v[74] = x[C.V74] * y[V.a_MTOR] / (x[C.K74] + y[V.a_MTOR])
        v[75] = x[C.V75] * y[V.a_AKT3_AKT1_AKT2] * y[V.BAD] / (x[C.K75] + y[V.BAD])
        v[76] = x[C.V76] * y[V.a_BAD] / (x[C.K76] + y[V.a_BAD])
        v[77] = x[C.V77] * y[V.a_AKT3_AKT1_AKT2] * y[V.GSK3B] / (x[C.K77] + y[V.GSK3B])
        v[78] = x[C.V78] * y[V.a_GSK3B] / (x[C.K78] + y[V.a_GSK3B])
        v[79] = x[C.V79] * y[V.a_AKT3_AKT1_AKT2] * y[V.CDKN1B] / (x[C.K79] + y[V.CDKN1B])
        v[80] = x[C.V80] * y[V.a_CDKN1B] / (x[C.K80] + y[V.a_CDKN1B])
        v[81] = x[C.V81] * y[V.a_AKT3_AKT1_AKT2] * y[V.CDKN1A] / (x[C.K81] + y[V.CDKN1A])
        v[82] = x[C.V82] * y[V.a_CDKN1A] / (x[C.K82] + y[V.a_CDKN1A])
        v[83] = x[C.V83] * y[V.a_NCK1_NCK2] * y[V.PAK4_BUB1B_PAK6_PAK1_PAK2_PAK3_PAK6_PAK5] / (x[C.K83] + y[V.PAK4_BUB1B_PAK6_PAK1_PAK2_PAK3_PAK6_PAK5])
        v[84] = x[C.V84] * y[V.a_PAK4_BUB1B_PAK6_PAK1_PAK2_PAK3_PAK6_PAK5] / (x[C.K84] + y[V.a_PAK4_BUB1B_PAK6_PAK1_PAK2_PAK3_PAK6_PAK5])
        v[85] = x[C.V85] * y[V.a_SOS1_SOS2] * y[V.HRAS_KRAS_NRAS] / (x[C.K85] + y[V.HRAS_KRAS_NRAS])
        v[86] = x[C.V86] * y[V.a_HRAS_KRAS_NRAS] / (x[C.K86] + y[V.a_HRAS_KRAS_NRAS])
        v[87] = x[C.V87] * y[V.a_CRK_CRKL] * y[V.ABL1_ABL2] / (x[C.K87] + y[V.ABL1_ABL2])
        v[88] = x[C.V88] * y[V.a_ABL1_ABL2] / (x[C.K88] + y[V.a_ABL1_ABL2])
        v[89] = x[C.V89] * y[V.a_SRC] * y[V.PTK2] / (x[C.K89] + y[V.PTK2])
        v[90] = x[C.V90] * y[V.a_PTK2] / (x[C.K90] + y[V.a_PTK2])
        v[91] = x[C.V91] * y[V.a_HRAS_KRAS_NRAS] * y[V.ARAF_RAF1_BRAF] / (x[C.K91] + y[V.ARAF_RAF1_BRAF])
        v[92] = x[C.V92] * y[V.a_ARAF_RAF1_BRAF] / (x[C.K92] + y[V.a_ARAF_RAF1_BRAF])
        v[93] = x[C.V93] * y[V.a_MTOR] * y[V.RPS6KB1_RPS6KB2] / (x[C.K93] + y[V.RPS6KB1_RPS6KB2])
        v[94] = x[C.V94] * y[V.a_RPS6KB1_RPS6KB2] / (x[C.K94] + y[V.a_RPS6KB1_RPS6KB2])
        v[95] = x[C.V95] * y[V.a_MTOR] * y[V.EIF4EBP1] / (x[C.K95] + y[V.EIF4EBP1])
        v[96] = x[C.V96] * y[V.a_EIF4EBP1] / (x[C.K96] + y[V.a_EIF4EBP1])
        v[97] = x[C.V97] * y[V.a_PAK4_BUB1B_PAK6_PAK1_PAK2_PAK3_PAK6_PAK5] * y[V.MAP2K7_MAP2K4] / (x[C.K97] + y[V.MAP2K7_MAP2K4])
        v[98] = x[C.V98] * y[V.a_MAP2K7_MAP2K4] / (x[C.K98] + y[V.a_MAP2K7_MAP2K4])
        v[99] = x[C.V99] * y[V.a_ARAF_RAF1_BRAF] * y[V.MAP2K1_MAP2K2] / (x[C.K99] + y[V.MAP2K1_MAP2K2])
        v[100] = x[C.V100] * y[V.a_MAP2K1_MAP2K2] / (x[C.K100] + y[V.a_MAP2K1_MAP2K2])
        v[101] = x[C.V101] * y[V.a_MAP2K7_MAP2K4] * y[V.MAPK8_MAPK9_MAPK10] / (x[C.K101] + y[V.MAPK8_MAPK9_MAPK10])
        v[102] = x[C.V102] * y[V.a_MAPK8_MAPK9_MAPK10] / (x[C.K102] + y[V.a_MAPK8_MAPK9_MAPK10])
        v[103] = x[C.V103] * y[V.a_MAP2K1_MAP2K2] * y[V.MAPK1_MAPK3] / (x[C.K103] + y[V.MAPK1_MAPK3])
        v[104] = x[C.V104] * y[V.a_MAPK1_MAPK3] / (x[C.K104] + y[V.a_MAPK1_MAPK3])
        v[105] = x[C.V105] * y[V.a_MAPK8_MAPK9_MAPK10] * y[V.JUN] / (x[C.K105] + y[V.JUN])
        v[106] = x[C.V106] * y[V.a_JUN] / (x[C.K106] + y[V.a_JUN])
        v[107] = x[C.V107] * y[V.a_MAPK8_MAPK9_MAPK10] * y[V.ELK1] / (x[C.K107] + y[V.ELK1])
        v[108] = x[C.V108] * y[V.a_ELK1] / (x[C.K108] + y[V.a_ELK1])
        v[109] = x[C.V109] * y[V.a_MAPK1_MAPK3] * y[V.MYC] / (x[C.K109] + y[V.MYC])
        v[110] = x[C.V110] * y[V.a_MYC] / (x[C.K110] + y[V.a_MYC])
        v[111] = x[C.V111] * y[V.a_MAPK1_MAPK3] * y[V.ELK1] / (x[C.K111] + y[V.ELK1])
        v[112] = x[C.kf112] * y[V.a_EGFR_EGFR]
        v[113] = x[C.kf113] * y[V.a_AKT3_AKT1_AKT2]
        v[114] = x[C.kf114] * y[V.a_MAPK1_MAPK3]
        v[115] = x[C.kf115] * y[V.a_MYC]

        return v
