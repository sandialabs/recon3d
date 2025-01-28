from pathlib import Path

RECON3D_PATH = Path(__file__).parent.parent.parent

BINARY_TO_SEMANTIC_YML = RECON3D_PATH.joinpath(
    "docs", "userguide", "src", "binary_to_semantic", "binary_to_semantic.yml"
)

DOWNSCALE_YML = RECON3D_PATH.joinpath(
    "docs", "userguide", "src", "downscale", "downscale.yml"
)

IMAGE_TO_VOXEL_YML = RECON3D_PATH.joinpath(
    "docs", "userguide", "src", "voxel_to_image", "image_to_voxel.yml"
)
INSTANCE_ANALYSIS_YML = RECON3D_PATH.joinpath(
    "docs", "userguide", "src", "instance_analysis", "instance_analysis.yml"
)

LETTER_F_3D_EXO = RECON3D_PATH.joinpath(
    "docs", "userguide", "src", "voxel_to_mesh", "letter_f_3d.exo"
)

LETTER_F_3D_YML = RECON3D_PATH.joinpath(
    "docs", "userguide", "src", "voxel_to_mesh", "letter_f_3d.yml"
)

SEMANTIC_TO_BINARY_YML = RECON3D_PATH.joinpath(
    "docs", "userguide", "src", "binary_to_semantic", "semantic_to_binary.yml"
)

VOXEL_TO_IMAGE_YML = RECON3D_PATH.joinpath(
    "docs", "userguide", "src", "voxel_to_image", "voxel_to_image.yml"
)
