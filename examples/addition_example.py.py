
from sasmodels.core import load_model_info
from sasmodels.sasview_model import make_model_from_info

model_info = load_model_info('fractal+mono_gauss_coil')
model_info.name = 'test2'
model_info.description = 'fractal + mono_gauss_coil'
Model = make_model_from_info(model_info)
