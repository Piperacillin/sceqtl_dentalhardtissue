#!/bin/bash

# ================= 全局配置中心 =================
# 在这里修改列表，00脚本和04脚本都会自动同步

EXPOSURE_LIST=(
    "onek1k"
    "1mscblood"
)

EXPOSURE_VAL_LIST=(
    "CD4_Dynamic"
    "dice"
)

OUTCOME_LIST=(
    "K11_ABRASION"
    "K11_EMBIMPACT_TEETH"
    "K11_EROSION_INCLAVO"
    "K11_ATTRITION"
    "K11_HYPO_ONLY"
    "K11_SUPNUM_ONLY_INCLAVO"
    "K11_CARIES_1_OPER_ONLYAVO"
    "K11_HYPOPLASIA_ENAMEL_INCLAVO"
    "K11_MIH"
    "K11_RESORPTION"
)

# 其他全局变量也可以放在这，比如
export GLOBAL_MAX_JOBS=69
