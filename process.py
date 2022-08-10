import datetime
from level1.write_lev1_nc import lev1_to_nc
from level2.write_lev2_nc import lev2_to_nc
from plots.generate_plots import generate_figure

start_date = datetime.date(2022, 7, 19)
end_date = datetime.date(2022, 7, 19)
delta = datetime.timedelta(days=1)
data_in = '/data/obs/site/jue/tophat/l1/'
# data_in = '/data/obs/site/cgn/foghat/l1/'
data_out_l1 = '/data/obs/site/jue/tophat/actris/level1/'
data_out_l2 = '/data/obs/site/jue/tophat/actris/level2/'

while start_date <= end_date:
    print(start_date)
    prev = start_date - delta
    nex = start_date + delta
    
    "Level 1"
    # lev1_to_nc('juelich', '1C01', data_in+start_date.strftime('%Y/%m/%d/'), data_in+prev.strftime('%Y/%m/%d/'), data_in+nex.strftime('%Y/%m/%d/'),
    #            data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc')
    
    "Level 2"
    # lev2_to_nc('juelich', '2I01', data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc', data_out_l2+start_date.strftime('%Y%m%d')+'_lwp.nc')
    lev2_to_nc('juelich', '2I02', data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc', data_out_l2+start_date.strftime('%Y%m%d')+'_iwv.nc')
    # lev2_to_nc('juelich', '2P02', data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc', data_out_l2+start_date.strftime('%Y%m%d')+'_tem_bl.nc')
    # lev2_to_nc('juelich', '2P03', data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc', data_out_l2+start_date.strftime('%Y%m%d')+'_hum.nc')
    # lev2_to_nc('juelich', '2P04', data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc', data_out_l2+start_date.strftime('%Y%m%d')+'_rh.nc')
    # lev2_to_nc('juelich', '2P07', data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc', data_out_l2+start_date.strftime('%Y%m%d')+'_pot.nc')    
    # lev2_to_nc('juelich', '2P08', data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc', data_out_l2+start_date.strftime('%Y%m%d')+'_eq_pot.nc')
    # lev2_to_nc('juelich', '2S02', data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc', data_out_l2+start_date.strftime('%Y%m%d')+'_tbx.nc')
    
    "Plot"
    # generate_figure(data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc',['tb'], ele_range = [89., 91.], save_path=data_out_l1, image_name='tb')
    # generate_figure(data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc',['ele', 'azi'], save_path=data_out_l1, image_name='sen')
    # generate_figure(data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc',['quality_flag'], save_path=data_out_l1, image_name='quality_flag')
    # generate_figure(data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc',['irt'], ele_range = [89., 91.], save_path=data_out_l1, image_name='irt')
    # generate_figure(data_out_l1+start_date.strftime('%Y%m%d')+'_lev1.nc', ['air_temperature','relative_humidity','air_pressure','rain_rate','wind_direction','wind_speed'], save_path=data_out_l1, image_name='met')

    # generate_figure(data_out_l2+start_date.strftime('%Y%m%d')+'_lwp.nc',['lwp'], save_path=data_out_l2, image_name='lwp')
    generate_figure(data_out_l2+start_date.strftime('%Y%m%d')+'_iwv.nc',['iwv'], save_path=data_out_l2, image_name='iwv')
    # generate_figure(data_out_l2+start_date.strftime('%Y%m%d')+'_tem_bl.nc',['temperature'], save_path=data_out_l2, image_name='tem_bl')
    # generate_figure(data_out_l2+start_date.strftime('%Y%m%d')+'_hum.nc',['water_vapor_vmr'], save_path=data_out_l2, image_name='hum')
    # generate_figure(data_out_l2+start_date.strftime('%Y%m%d')+'_rh.nc',['relative_humidity'], save_path=data_out_l2, image_name='rh')
    # generate_figure(data_out_l2+start_date.strftime('%Y%m%d')+'_pot.nc',['potential_temperature'], save_path=data_out_l2, image_name='pot')
    # generate_figure(data_out_l2+start_date.strftime('%Y%m%d')+'_eq_pot.nc',['equivalent_potential_temperature'], ele_range = [89., 91.], save_path=data_out_l2, image_name='eq_pot')
    # generate_figure(data_out_l2+start_date.strftime('%Y%m%d')+'_tbx.nc',['tb_spectrum'], save_path=data_out_l2, image_name='tbx')
    
    start_date += delta
