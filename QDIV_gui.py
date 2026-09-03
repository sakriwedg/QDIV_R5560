# -*- coding: utf-8 -*-
"""
QDIV-gui : Python Graphical User Interface for Position-Sensitive Neutron Detectors read-out by the resistive-charge-division method implemented in the QDIV firmware initially developped at the ILL and running in CAEN digitizers

Created on Sat Jul 11 07:41:13 2026

@author: marchalj
"""
version="23.07.2026"

import plotly.express as px
import asyncio
from nicegui import app,run,ui,html
import numpy as np
import time
from datetime import date
#import pandas as pd
from plotly.subplots import make_subplots
import csv
import math
import os
import random
import threading

######################################################################
# GUI iniitialisation
######################################################################
# QDIV-gui modules
import detectors
import config
import caen

# Global variables
RMS_TABLE=[]
OFFSET_TABLE=[]
TRACES_COUNTER=0
TRACES_TO_SAVE_COUNTER=0
global IM,PHS,COUNTS,LIST_DATA
IM=np.zeros([config.det.nb_of_tubes,config.nb_of_pixels_per_tube]) # Position histogram (2D image)
PHS=np.zeros([config.det.nb_of_tubes,config.nb_of_bins_in_spectrum]) # A+B sum histogram (Pulse Height Spectra)
PH_VS_POSITION=np.zeros([config.det.nb_of_tubes,config.nb_of_pixels_per_tube]) # Pulse Height vs position 
PHS_A=PHS
PHS_B=PHS
COUNTS=0
LIST_DATA=[]
ELAPSED_TIME=1
# Data read-out thread : the digitizer buffers are emptied in a dedicated thread so that the
# plot refresh (which runs in the GUI thread) never interrupts the read-out and the count rate
# is preserved
DATA_LOCK=threading.Lock() # protects the histograms shared by the read-out thread and the plots
DATA_READER_THREAD=None
DATA_READER_STOP_EVENT=None
DATA_READER_ERROR=None
LAST_BUFFER_READ_TIME=0
LIST_SAVING_ENABLED=False # snapshot of list_saving_checkbox (the thread must not touch the GUI)

# Virtual digitizer class for GUI debug mode
class debug_mode_digitizer:
    def __init__(self):            
        self.ADC_sampling_freq=120E6 #MHz
        self.ADC_input_range=2 #V
        self.ADC_bit_depth=14
        self.filter_output_bit_depth=16
        self.mV_per_ADU=self.ADC_input_range*1000/2**self.ADC_bit_depth
        self.custom_packet_buffer_size=1024
    
    def set_registers(self,registers_table):
        pass

    def reset_pulse_height_spectrum(self,tube_number):
        pass

    def read_pulse_height_spectrum(self,tube_number):
        bins = []
        phs = []
        buf_info_valid_bins=self.custom_packet_buffer_size
        buf_data=np.random.randint(100,size=buf_info_valid_bins)
        for index in range(buf_info_valid_bins):
            bins.append(index)
            phs.append(buf_data[index])
        return(bins,np.array(phs))

    def read_scope_traces(self,scope_settings):
        channels=4
        samples=self.custom_packet_buffer_size
        analog_traces = np.zeros([channels,samples])
        digital_traces = np.zeros([channels,samples])
        for i in range(4):
            analog_traces[i]=np.random.poisson(lam=0.1+i,size=samples)+i
            digital_traces[i]=np.random.randint(2,size=samples)+i
        analog_time_axis=1E6*np.arange(self.custom_packet_buffer_size)/(self.ADC_sampling_freq)
        res_osc=0
        return(res_osc,analog_time_axis,analog_traces,digital_traces)
        
    def configure_data_read_out(self):
        pass
    
    def start_data_read_out(self):
        return time.time()

    def stop_data_read_out(self):
        return time.time()

    def read_data(self):
        nb_of_events=100
        channel_list=[]
        energy_A_list=[]
        energy_B_list=[]
        [channel_list.append(random.randint(0, config.det.nb_of_tubes-1)) for n in range(nb_of_events)]
        [energy_A_list.append(random.randint(0, 2**self.filter_output_bit_depth-1)) for n in range(nb_of_events)]
        [energy_B_list.append(random.randint(0, 2**self.filter_output_bit_depth-1)) for n in range(nb_of_events)]
        buffer_read_time=time.time()
        return(channel_list,energy_A_list,energy_B_list,buffer_read_time)
    
# Digitizer object creation
if config.debug_mode:
    digitizer=debug_mode_digitizer()
else:
    digitizer=caen.digitizer()
    
# Data saving folder creation 
def check_data_folder(path):
        if os.path.exists(path) == 0:
            os.makedirs(path)
            print('Directory ' + path + ' created')

######################################################################
# Data acquisition
######################################################################

def reset_acquisition_variables():
    global IM,PHS,PHS_A,PHS_B,COUNTS,LIST_DATA,PH_VS_POSITION
    with DATA_LOCK:
        IM=np.zeros([config.det.nb_of_tubes,config.nb_of_pixels_per_tube])
        PHS=np.zeros([config.det.nb_of_tubes,config.nb_of_bins_in_spectrum])
        PHS_A=PHS
        PHS_B=PHS
        COUNTS=0
        LIST_DATA=[]
        PH_VS_POSITION=np.zeros([config.det.nb_of_tubes,config.nb_of_pixels_per_tube])

def read_data():
    global IM,PHS,PHS_A,PHS_B,COUNTS,LIST_DATA,PH_VS_POSITION
    [channel_list,energy_A_list,energy_B_list,buffer_read_time]=digitizer.read_data()


    # check if energy_A_list and energy_B_list are not empty to avoid division by zero
    #if len(energy_A_list) == 0 or len(energy_B_list) == 0:
    #    print("No data read from digitizer.")
    # check if channel_list is not empty to avoid errors in histogramming
    #if len(channel_list) == 0:
    #    print("No channel data read from digitizer.")
    #else:
    #print(f"Read {len(channel_list)} events: {len(energy_A_list)} energy A values, {len(energy_B_list)} energy B values.")

    # check whether energy_A + energy_B is zero to avoid division by zero in position calculation
    if any((np.array(energy_A_list) + np.array(energy_B_list)) == 0):

        # count number of events with zero total energy
        zero_energy_events = np.sum((np.array(energy_A_list) + np.array(energy_B_list)) == 0)
        #print(f"Warning: {zero_energy_events} events have zero total energy, which may lead to division by zero in position calculation.")
        #print("Warning: Some events have zero total energy, which may lead to division by zero in position calculation.")
        # Optionally, you can filter out these events or handle them appropriately
        valid_indices = (np.array(energy_A_list) + np.array(energy_B_list)) != 0
        channel_list = np.array(channel_list)[valid_indices]
        energy_A_list = np.array(energy_A_list)[valid_indices]
        energy_B_list = np.array(energy_B_list)[valid_indices]


    channel_array=np.array(channel_list,dtype=np.int32)
    energy_A_array=np.array(energy_A_list,dtype=np.float64)
    energy_B_array=np.array(energy_B_list,dtype=np.float64)
    total_energy_array=energy_A_array+energy_B_array
    # Reject invalid events : channel number outside of the detector (garbage words still
    # sitting in the digitizer buffers) or null total energy (division by zero in the
    # charge-division position calculation)
    valid_array=(channel_array>=0)&(channel_array<config.det.nb_of_tubes)&(total_energy_array>0)
    nb_of_rejected_events=int(np.size(valid_array)-np.count_nonzero(valid_array))
    if nb_of_rejected_events and config.chatty:
        print('######## Rejected events : '+str(nb_of_rejected_events))
    channel_array=channel_array[valid_array]
    energy_A_array=energy_A_array[valid_array]
    energy_B_array=energy_B_array[valid_array]
    total_energy_array=total_energy_array[valid_array]
    # Charge division : position is in [0,nb_of_pixels_per_tube-1] (B=0 would give an extra bin)
    position_array=np.minimum(np.int32(config.nb_of_pixels_per_tube*energy_A_array/total_energy_array),config.nb_of_pixels_per_tube-1)
    yrange=(np.arange(0,2*2**16+1,int(2*2**16/config.nb_of_bins_in_spectrum)))-0.5
    xrange=np.arange(config.det.nb_of_tubes+1)-0.5
    HISTO_ENERGY=np.histogram2d(channel_array,total_energy_array,bins=[xrange,yrange])
    HISTO_ENERGY_A=np.histogram2d(channel_array,energy_A_array,bins=[xrange,yrange])
    HISTO_ENERGY_B=np.histogram2d(channel_array,energy_B_array,bins=[xrange,yrange])
    yrange=np.arange(config.nb_of_pixels_per_tube+1)-0.5
    HISTO2D=np.histogram2d(channel_array,position_array,bins=[xrange,yrange])
    HISTO2D[0],HISTO_ENERGY[0],HISTO_ENERGY_A[0],HISTO_ENERGY_B[0]
    energy_vs_2D_pos_array=np.zeros([config.det.nb_of_tubes,config.nb_of_pixels_per_tube])
    counts_vs_2D_pos_array=np.zeros([config.det.nb_of_tubes,config.nb_of_pixels_per_tube])
    np.add.at(energy_vs_2D_pos_array,(channel_array,position_array),total_energy_array)
    np.add.at(counts_vs_2D_pos_array,(channel_array,position_array),1)
    PH_vs_position=np.divide(energy_vs_2D_pos_array,counts_vs_2D_pos_array,out=np.zeros_like(energy_vs_2D_pos_array),where=counts_vs_2D_pos_array>0)
    # Accumulation of the histograms shared with the plots (this function runs in the read-out
    # thread : everything above works on local arrays only, the lock is taken here and released
    # immediately so that the read-out is never slowed down by a plot refresh)
    with DATA_LOCK:
        COUNTS=COUNTS+np.sum(HISTO2D[0])
        IM=IM+HISTO2D[0]
        PHS=PHS+HISTO_ENERGY[0]
        PHS_A=PHS_A+HISTO_ENERGY_A[0]
        PHS_B=PHS_B+HISTO_ENERGY_B[0]
        PH_VS_POSITION=PH_VS_POSITION+PH_vs_position
        if LIST_SAVING_ENABLED:
            # Create a 2D array with 3 columns
            combined_array = np.column_stack((channel_array,energy_A_array ,energy_B_array))
            # Convert to list and append multiple times
            LIST_DATA.extend(combined_array.tolist())
    return buffer_read_time,int(np.size(valid_array))

def data_reader_loop(stop_event):
    """Empty the digitizer buffers as fast as possible until stop_event is set.

    Runs in its own thread : the GIL is released during the SciSDK read and during the numpy
    histogramming, so the plot refresh of the GUI thread no longer stops the data flow.
    """
    global LAST_BUFFER_READ_TIME,DATA_READER_ERROR
    nb_of_empty_reads=0
    while not stop_event.is_set():
        try:
            buffer_read_time,nb_of_events=read_data()
        except Exception as error:
            DATA_READER_ERROR=error
            print('!!! Data read-out thread stopped : '+repr(error))
            return
        LAST_BUFFER_READ_TIME=buffer_read_time
        if nb_of_events:
            nb_of_empty_reads=0
        else:
            nb_of_empty_reads=nb_of_empty_reads+1
            if nb_of_empty_reads>10:
                # Read-out is non-blocking : when nothing comes any more, back off instead of
                # polling flat out, otherwise the polling eats a full CPU core and slows the GUI
                # down. A single empty buffer is never waited on, so no count rate is lost.
                stop_event.wait(10E-3)

def start_data_reader_thread():
    global DATA_READER_THREAD,DATA_READER_STOP_EVENT,DATA_READER_ERROR,LAST_BUFFER_READ_TIME
    stop_data_reader_thread()
    DATA_READER_ERROR=None
    LAST_BUFFER_READ_TIME=time.time()
    DATA_READER_STOP_EVENT=threading.Event()
    DATA_READER_THREAD=threading.Thread(target=data_reader_loop,args=(DATA_READER_STOP_EVENT,),name='QDIV data read-out',daemon=True)
    DATA_READER_THREAD.start()

def stop_data_reader_thread():
    global DATA_READER_THREAD,DATA_READER_STOP_EVENT
    if DATA_READER_THREAD is not None:
        DATA_READER_STOP_EVENT.set()
        DATA_READER_THREAD.join(timeout=5)
        if DATA_READER_THREAD.is_alive():
            print('!!! Data read-out thread did not stop within 5s')
        DATA_READER_THREAD=None

def purge_buffers(nb_of_reads=20):
    [read_data() for n in range(nb_of_reads)]

def update_acquisition_figures():
    global IM,PHS,PHS_A,PHS_B,COUNTS,LIST_DATA,ELAPSED_TIME,PH_VS_POSITION
    # Snapshot of the histograms being filled by the read-out thread : the figures are then
    # built outside of the lock, so plotting never blocks the read-out
    with DATA_LOCK:
        IM_snapshot=IM.copy()
        PHS_snapshot=PHS.copy()
        PHS_A_snapshot=PHS_A.copy()
        PHS_B_snapshot=PHS_B.copy()
        PH_VS_POSITION_snapshot=PH_VS_POSITION.copy()
        COUNTS_snapshot=COUNTS

    match log_checkbox.value:
        case False:
            IM_fig=IM_snapshot
            PHS_fig=PHS_snapshot
            PH_vs_position_fig=PH_VS_POSITION_snapshot
        case True:
            IM_fig=np.log(IM_snapshot)
            PHS_fig=np.log(PHS_snapshot)
            PH_vs_position_fig=np.log(PH_VS_POSITION_snapshot)

    counts_per_second=str(int(COUNTS_snapshot/ELAPSED_TIME))+ ' counts/s'

    match acq_fig_selec_toggle.value:
        case 'Image':
            ACQ_FIG1.figure=px.imshow(IM_fig,aspect=config.det.aspect_ratio,labels=dict(x="X channel", y="Y channel", color="Counts"))
            ACQ_FIG1.figure.update_layout(title={'text':counts_per_second,'x': 0.5})
            ACQ_FIG1.figure.update_xaxes(title='Tube length (a.u.)')
            ACQ_FIG1.figure.update_yaxes(title='Tube number')
            ACQ_FIG1.set_visibility(True)
            ACQ_FIG1.update()
        case 'Pulse Height Spectra':
            match projection_checkbox.value:
                case False:
                    #print("NO PROJ")
                    ACQ_PROJ.set_visibility(False)
                    ACQ_FIG1.set_visibility(True)
                    ACQ_FIG1.figure=px.imshow(PHS_fig,aspect=config.det.aspect_ratio,labels=dict(x="X channel", y="Y channel", color="Counts"))
                    ACQ_FIG1.figure.update_layout(title={'text':counts_per_second,'x': 0.5})
                    ACQ_FIG1.figure.update_xaxes(title='Bins')
                    ACQ_FIG1.figure.update_yaxes(title='Tube number')
                    #ACQ_FIG1.classes('w-full justify-center no-wrap')
                case True:
                    print(" PROJ")
                    ACQ_FIG1.set_visibility(False)
                    ACQ_PROJ.set_visibility(True)
                    ax = ACQ_PROJ.figure.gca()
                    for i, row in enumerate(PHS_snapshot):
                        ax.plot(row)
                    ax.tick_params(axis='x', labelsize=6) 
                    ax.tick_params(axis='y', labelsize=6)
                    ax.set_xlabel("Bins",fontsize=5)
                    ax.set_ylabel("Counts", fontsize=5)
                    #if acq_fig_selec_toggle.value:
                    #     ax.update(yaxis_type="log")
            ACQ_PROJ.update()
            ACQ_FIG1.update()
        case 'Single tube':
            ACQ_FIG2.push(range(config.nb_of_pixels_per_tube),[IM_snapshot[int(tube_number.value-1),:]])
            ACQ_FIG1.set_visibility(False)
            ACQ_FIG3.push(range(config.nb_of_bins_in_spectrum),[PHS_snapshot[int(tube_number.value-1),:],PHS_A_snapshot[int(tube_number.value-1),:],PHS_B_snapshot[int(tube_number.value-1),:]])
            # if acq_fig_selec_toggle.value:
            #     ACQ_FIG2.fig.update(yaxis_type="log")
            #     ACQ_FIG3.fig.update(yaxis_type="log")
            #TODO : log for 1D plots
            ACQ_FIG3.update()
            ACQ_FIG2.update()
            ACQ_FIG1.update()
        case 'Pulse height vs position':
            ACQ_FIG1.figure=px.imshow(PH_vs_position_fig,aspect=config.det.aspect_ratio,labels=dict(x="X channel", y="Y channel", color="Counts"))
            ACQ_FIG1.figure.update_layout(title={'text':counts_per_second,'x': 0.5})
            ACQ_FIG1.figure.update_xaxes(title='Tube length [a.u.]')
            ACQ_FIG1.figure.update_yaxes(title='Counts')
            ACQ_FIG1.update()

def save_acquisition_files():
            fullPath=path.value+'/'+subpath.value
            if os.path.isfile(fullPath+'/'+'idx.txt'):
                fileIndex = np.loadtxt(fullPath+'/'+'idx.txt', dtype='int16')
                fileIndex = fileIndex+1
            else:
                fileIndex = 1
                np.savetxt(fullPath+'/'+'idx.txt',np.array([fileIndex], np.int32), fmt="%05d")

            dataFile=fullPath+'/'+'image'+'_'+str("%05d" % fileIndex)+'.txt'
            np.savetxt(dataFile, IM , fmt='%i')
            dataFile=fullPath+'/'+'PHS'+'_'+str("%05d" % fileIndex)+'.txt'
            np.savetxt(dataFile, PHS , fmt='%i')
            dataFile=fullPath+'/'+'PH_vs_position'+'_'+str("%05d" % fileIndex)+'.txt'
            np.savetxt(dataFile, PH_VS_POSITION , fmt='%i')

            acq_fig_selec_toggle.value='Image'
            update_acquisition_figures()
            dataFile=fullPath+'/'+'image'+'_'+str("%05d" % fileIndex)+'.pdf'
            ACQ_FIG1.figure.write_image(dataFile)  
            acq_fig_selec_toggle.value='Pulse Height Spectra'
            update_acquisition_figures()
            dataFile=fullPath+'/'+'PHS'+'_'+str("%05d" % fileIndex)+'.pdf'
            ACQ_FIG1.figure.write_image(dataFile) 
            acq_fig_selec_toggle.value='Pulse height vs position'
            update_acquisition_figures()
            dataFile=fullPath+'/'+'PH_vs_position'+'_'+str("%05d" % fileIndex)+'.pdf'
            ACQ_FIG1.figure.write_image(dataFile) 

            if list_saving_checkbox.value:
                dataFile=fullPath+'/'+'list_data'+'_'+str("%05d" % fileIndex)+'.txt'
                np.savetxt(dataFile, np.array(LIST_DATA) , fmt='%i')

            ui.notification("Plots saved in .txt and .pdf formats",timeout=1)

            if config.chatty:
                print('### Traces saved in : ' + dataFile)

async def start_acquisition():
    global IM,PHS,PHS_A,PHS_B,COUNTS,LIST_DATA,ELAPSED_TIME,GAIN_MAP,LIST_SAVING_ENABLED
    print("Acquisition started")
    ui.update()

    registers_table_update()
    digitizer.set_registers(table)
    digitizer.start_data_read_out()
    digitizer.stop_data_read_out()
    digitizer.configure_data_read_out()
    reset_acquisition_variables()
    LIST_SAVING_ENABLED=list_saving_checkbox.value # the read-out thread must not read the GUI
    await run.io_bound(purge_buffers) # Dummy read-outs to purge buffers
    notif = ui.notification(timeout=None,spinner=True,position="bottom-left")
    nIm=0
    if startToggleButton._state:
        while nIm < (np.int16(repetitions.value)) and startToggleButton._state==1 and DATA_READER_ERROR is None:
            frameNb = 1
            reset_acquisition_variables()
            LIST_SAVING_ENABLED=list_saving_checkbox.value
            start_time=digitizer.start_data_read_out()
            # The buffers are emptied by a dedicated thread : the plot refresh below runs in
            # parallel with the read-out instead of stopping it (the count rate was lost during
            # every figure update)
            start_data_reader_thread()
            while time.time() < start_time  + expT.value and startToggleButton._state==1:
                if startToggleButton._state:
                    await asyncio.sleep(max(0,start_time+frameNb*(config.GUI_refresh_period)-time.time()))
                    ELAPSED_TIME=max(LAST_BUFFER_READ_TIME-start_time,1E-3)
                    update_acquisition_figures()
                    if not live_mode_switch.value:
                        notif.message = 'Acquisition : '+ str(nIm+1)+ '/'+str(np.int16(repetitions.value)) + ' ; Exposure time : ' + str(np.int16(time.time()-start_time ))+'s/'+str(np.int16(expT.value))+'s'
                    frameNb = frameNb+1
                    if DATA_READER_ERROR is not None:
                        ui.notify('Data read-out stopped : '+repr(DATA_READER_ERROR),type='negative')
                        break
            await run.io_bound(stop_data_reader_thread)
            stop_time=digitizer.stop_data_read_out()
            await run.io_bound(purge_buffers)
            ELAPSED_TIME=stop_time-start_time
            update_acquisition_figures()
            nIm=nIm+1
            if plots_saving_checkbox.value and startToggleButton._state:
                check_data_folder(path.value+subpath.value )
                save_acquisition_files()

    else:
        pass
    stop_data_reader_thread()
    notif.dismiss()
    ui.update
    return

async def live_mode():
    while live_mode_switch.value:
        acq_Trigger_Thresh.set_visibility(False)
        startToggleButton.set_visibility(False)
        repetitions.set_visibility(False)
        plots_saving_checkbox.set_visibility(False)
        ui.update()
        await startToggleButton.toggle()
    startToggleButton.set_visibility(True)
    repetitions.set_visibility(True)
    plots_saving_checkbox.set_visibility(True)
    acq_Trigger_Thresh.set_visibility(True)

######################################################################
# Oscilloscope
######################################################################

def push_traces_line_plots(time_axis,analog_traces,digital_traces):
    global RMS_TABLE, OFFSET_TABLE
    analog_input_line_plot.push(time_axis, [analog_traces[0]*digitizer.mV_per_ADU-digitizer.ADC_input_range*1000/2, analog_traces[1]*digitizer.mV_per_ADU-digitizer.ADC_input_range*1000/2])
    analog_input_line_plot.fig.suptitle("Filter input signals (Tube "+str(tube_number.value)+ ')')
    analog_output_line_plot.push(time_axis, [analog_traces[2], analog_traces[3]])
    analog_output_line_plot.fig.suptitle("Filter output signals (Tube "+str(tube_number.value)+ ')')
    
    if stats_or_spectra_toggle.value=='Digital signals':
        digital_line_plot.push(time_axis, [digital_traces[0], digital_traces[1],digital_traces[2], digital_traces[3]])

    if stats_or_spectra_toggle.value=='Statistics':
        rms_val_ch1.value=(round(np.std(analog_traces[0]),2))*digitizer.mV_per_ADU
        rms_val_ch2.value=(round(np.std(analog_traces[1]),2))*digitizer.mV_per_ADU
        rms_val_ch3.value=(round(np.std(analog_traces[2]),2))
        rms_val_ch4.value=(round(np.std(analog_traces[3]),2))

        mean_val_ch1.value=np.mean(analog_traces[0])*digitizer.mV_per_ADU-digitizer.ADC_input_range*1000/2
        mean_val_ch2.value=np.mean(analog_traces[1])*digitizer.mV_per_ADU-digitizer.ADC_input_range*1000/2
        mean_val_ch3.value=np.mean(analog_traces[2])
        mean_val_ch4.value=np.mean(analog_traces[3])

    RMS_TABLE.append([rms_val_ch1.value,rms_val_ch2.value,rms_val_ch3.value,rms_val_ch4.value])
    OFFSET_TABLE.append([mean_val_ch1.value,mean_val_ch2.value,mean_val_ch3.value,mean_val_ch4.value])

def save_traces(combined_array):
        fullPath=traces_path.value+'/'+traces_subpath.value 
        if os.path.isfile(fullPath+'idx.txt'):
            fileIndex = np.loadtxt(fullPath+'idx.txt', dtype='int16')
            fileIndex = fileIndex+1
        else:
            fileIndex = 1
        np.savetxt(fullPath+'idx.txt',np.array([fileIndex], np.int32), fmt="%05d")
        dataFile=fullPath+traces_filename.value+'tube#'+str(round(tube_number.value-1))+'_'+str("%05d" % fileIndex)+'.txt'
        np.savetxt(dataFile, combined_array  , fmt='%i')
        if config.chatty:
            print('### Traces saved in : ' + dataFile)

def update_scope_traces():
    global TRACES_TO_SAVE_COUNTER,TRACES_COUNTER,RMS_TABLE, OFFSET_TABLE
    check_data_folder(traces_path.value+'/'+traces_subpath.value )
    registers_table_update()
    digitizer.set_registers(table)
    scope_settings={'tube_number':tube_number.value,'threshold_slider':scope_threshold_slider.value,'trigger_mode':scope_trigger_mode.value,'trig_ch': scope_trig_ch_button.options.index(scope_trig_ch_button.value),'decimation_factor':decimation_factor.value}
    #buf_osc=digitizer.updateScopeSettings(sdk_list,scope_settings)
    [res_osc,analog_time_axis,analog_traces,digital_traces]=digitizer.read_scope_traces(scope_settings)
    if res_osc == 0:
        push_traces_line_plots(analog_time_axis,analog_traces,digital_traces)
    else:
        ui.notify('No trigger (check scope threshold)')
        ui.update()
        return
    if traces_autosave_checkbox.value:
        TRACES_TO_SAVE_COUNTER=TRACES_TO_SAVE_COUNTER+1 
        if TRACES_TO_SAVE_COUNTER<number_of_traces_input.value+1:
            combined_array = np.column_stack((analog_traces[0], analog_traces[1],analog_traces[2], analog_traces[3]))
            save_traces(combined_array)
            ui.notify("Saving traces", timeout=0.1)
        else:
            traces_autosave_checkbox.value=False
            scope_switch.value=False
            TRACES_TO_SAVE_COUNTER=0
        ui.update()
    TRACES_COUNTER=TRACES_COUNTER+1
    
    if stats_or_spectra_toggle.value=='A+B pulse height spectrum':
            (bins,phs)=digitizer.read_pulse_height_spectrum(tube_number.value)
            spectra_plot.push(bins, [phs])

######################################################################
# Electronic noise measurement
######################################################################

def fig_noise_and_offset_layout_vin(fig_noise_and_offset):
    fig_noise_and_offset[0].update_yaxes(title_text='mV rms')
    fig_noise_and_offset[1].update_yaxes(title_text='mV')
    fig_noise_and_offset[0].update_layout(title_text='Filter input', title_x=0.5)
    fig_noise_and_offset[1].update_layout(title_text='Filter input', title_x=0.5)
    for f in range(2):
        fig_noise_and_offset[f].update_xaxes(title_text='Tube number')
        series_names = ["Channel A", "Channel B"]
        for i, name in enumerate(series_names):
            fig_noise_and_offset[f].data[i].name = name
        fig_noise_and_offset[f].update_xaxes(
            tickmode="linear",   # graduations régulières
            tick0=0,              # première graduation
            dtick=1,              # pas de 1 unité
            tickformat="d"       # affichage sans décimales
        )
        # ---- Placement de la légende à l’intérieur ----
        fig_noise_and_offset[f].update_layout(legend_title_text='')
        fig_noise_and_offset[f].update_layout(
            legend=dict(
                # Position (0 = bord gauche/bas, 1 = bord droit/haut)
                x=0.98,          # légèrement à l’intérieur du bord gauche
                y=0.98,          # près du bord supérieur
                xanchor='right', # ancre horizontale
                yanchor='top',  # ancre verticale
                bgcolor='rgba(255,255,255,0.8)',  # fond semi‑transparent
                bordercolor='gray',
                borderwidth=1
            )
        )

def fig_noise_and_offset_layout_vout(fig_noise_and_offset):
    fig_noise_and_offset[0].update_yaxes(title_text='ADU rms')
    fig_noise_and_offset[1].update_yaxes(title_text='ADUs')
    fig_noise_and_offset[0].update_layout(title_text='Filter output', title_x=0.5)
    fig_noise_and_offset[1].update_layout(title_text='Filter output', title_x=0.5)
    for f in range(2):
        fig_noise_and_offset[f].update_xaxes(title_text='Tube number')
        series_names = ["Channel A", "Channel B"]
        for i, name in enumerate(series_names):
            fig_noise_and_offset[f].data[i].name = name
        fig_noise_and_offset[f].update_xaxes(
            tickmode="linear",   # graduations régulières
            tick0=0,              # première graduation
            dtick=1,              # pas de 1 unité
            tickformat="d"       # affichage sans décimales
        )
        # ---- Placement de la légende à l’intérieur ----
        fig_noise_and_offset[f].update_layout(legend_title_text='')
        fig_noise_and_offset[f].update_layout(
            legend=dict(
                # Position (0 = bord gauche/bas, 1 = bord droit/haut)
                x=0.98,          # légèrement à l’intérieur du bord gauche
                y=0.98,          # près du bord supérieur
                xanchor='right', # ancre horizontale
                yanchor='top',  # ancre verticale
                bgcolor='rgba(255,255,255,0.8)',  # fond semi‑transparent
                bordercolor='gray',
                borderwidth=1
            )
        )

def update_fig_noise_and_offset(fig_noise_vin,fig_noise_vout,fig_offset_vin,fig_offset_vout):
    match noise_or_offset_toggle.value:
        case 'RMS noise':
            FIG_NOISE_AND_OFFSET_VIN.figure=fig_noise_vin
            FIG_NOISE_AND_OFFSET_VIN.set_visibility(True)
            FIG_NOISE_AND_OFFSET_VIN.update()
            FIG_NOISE_AND_OFFSET_VOUT.figure=fig_noise_vout
            FIG_NOISE_AND_OFFSET_VOUT.set_visibility(True)
            FIG_NOISE_AND_OFFSET_VOUT.update()
        case 'Offsets':
            FIG_NOISE_AND_OFFSET_VIN.figure=fig_offset_vin
            FIG_NOISE_AND_OFFSET_VIN.set_visibility(True)
            FIG_NOISE_AND_OFFSET_VIN.update()
            FIG_NOISE_AND_OFFSET_VOUT.figure=fig_offset_vout
            FIG_NOISE_AND_OFFSET_VOUT.set_visibility(True)
            FIG_NOISE_AND_OFFSET_VOUT.update()
#TODO start tube numbering at 1

async def start_noise_measurement():
    global RMS_TABLE
    global OFFSET_TABLE
    global TRACES_COUNTER
    TRACES_COUNTER=0
    tube_RMS_mean=[]
    tube_OFFSET_mean=[]
    traces_subpath_copy=traces_subpath.value
    traces_subpath.value=noise_subpath.value
    scope_trigger_mode.value='self'
    check_data_folder(noise_path.value+'/'+noise_subpath.value )
    for tube in range(config.det.nb_of_tubes):
        tube_number.value=tube+1
        print("Tube number:", tube_number.value)
        stats_or_spectra_toggle.value='Statistics'
        if noise_meas_button._state:
            RMS_TABLE=[]
            OFFSET_TABLE=[]
            n=ui.notification("Saving oscilloscope traces for tube #" + str(tube) )
            traces_autosave_checkbox.value=True
            scope_switch.value=True
            while scope_switch.value==True:
                await asyncio.sleep(1)
            n.dismiss()
            scope_switch.value=False
            tube_RMS_mean.append(np.sqrt(np.sum(np.square(np.array(RMS_TABLE)),0)/np.shape(np.array(RMS_TABLE))[0]))
            tube_OFFSET_mean.append(np.mean(np.array(OFFSET_TABLE),0))
            fig_noise_vin=px.line(np.array(tube_RMS_mean)[:,:2],title='RMS noise per tube',markers=True)
            fig_offset_vin=px.line(np.array(tube_OFFSET_mean)[:,:2],title='Offset per tube',markers=True)
            fig_noise_and_offset_layout_vin([fig_noise_vin,fig_offset_vin])
            fig_noise_vout=px.line(np.array(tube_RMS_mean)[:,2:],title='RMS noise per tube VOUT',markers=True)
            fig_offset_vout=px.line(np.array(tube_OFFSET_mean)[:,2:],title='Offset per tube',markers=True)
            fig_noise_and_offset_layout_vout([fig_noise_vout,fig_offset_vout])
            update_fig_noise_and_offset(fig_noise_vin,fig_noise_vout,fig_offset_vin,fig_offset_vout)
            await asyncio.sleep(0.1)
        else:
            pass
    traces_subpath.value=traces_subpath_copy
    traces_autosave_checkbox.value=False
    fullPath=noise_path.value+'/'+noise_subpath.value 
    dataFile=fullPath+'rms_noise.txt'
    print("dataFile:", dataFile)
    np.savetxt(dataFile, tube_RMS_mean , fmt='%.2f')
    dataFile=fullPath+'offsets.txt'
    np.savetxt(dataFile, tube_OFFSET_mean , fmt='%.2f')
    dataFile=fullPath+'Vin_rms_noise.pdf'
    fig_noise_vin.write_image(dataFile)
    dataFile=fullPath+'Vin_offsets.pdf'
    fig_offset_vin.write_image(dataFile)
    dataFile=fullPath+'Vout_rms_noise.pdf'
    fig_noise_vout.write_image(dataFile)
    dataFile=fullPath+'Vout_offsets.pdf'
    fig_offset_vout.write_image(dataFile) 
    ui.update

######################################################################
# GUI functions
######################################################################

def update_GUI_scope_settings(): 

    def get_value(rows, key_search, value_search, key_return):
        """
        Retourne la valeur associée à `key_return` dans la première ligne où
        `key_search` == `value_search`. Retourne None si aucune ligne ne correspond.
        """
        for row in rows:
            if row.get(key_search) == value_search:
                #print("row.get(key_search):", row.get(key_search))
                #print("row['type']:", row['type'])
                match row['type']:
                    case 'Decimal':
                        #print("row.get(key_return):", row.get(key_return))
                        return int(row.get(key_return))
                    case 'Hexadecimal':
                        #print("int(row.get(key_return),16):", int(row.get(key_return),16))
                        return(int(row.get(key_return),16))
        return None

    if config.chatty:  
        print("Update signal processing parameters in the oscilloscope menu of the GUI using the registers values present in the register menu of the GUI")
    Trigger_Thresh.value=get_value(table.rows,'name', 'PCFG_THRS','value')
    #Trigger_Thresh.value=mV_per_ADU*get_value(table.rows,'name', 'PCFG_THRS','value')
    Trigger_Source.value=(Trigger_Source.options[(get_value(table.rows,'name', 'PCFG_TRG_SOURCE','value'))])
    Trigger_Mode.value=(Trigger_Mode.options[(get_value(table.rows,'name', 'PCFG_TRG_MODE','value'))])
    Trigger_hysteresis.value=get_value(table.rows,'name', 'PCFG_HIST','value')
    Trigger_hold_off.value=(1E6/digitizer.ADC_sampling_freq)*get_value(table.rows,'name', 'PCFG_TRG_HOLD','value')
    Trigger_output_width.value=(1E6/digitizer.ADC_sampling_freq)*get_value(table.rows,'name', 'PCFG_TRG_W','value')
    Gate_window.value =(1E6/digitizer.ADC_sampling_freq)*get_value(table.rows,'name', 'PCFG_GATE_W','value')
    polarity.value=(polarity.options[int(get_value(table.rows,'name','PCFG_POL','value'))])
    #offset.value=get_value(table.rows,'name', 'PCFG_RAWOFS','value')
    offset.value=digitizer.mV_per_ADU*get_value(table.rows,'name', 'PCFG_RAWOFS','value')
    # OFFSETS MENU TO BE MODIFIED SINCE OFFSETS SEEM TOP BE PER CHANNEL IN THIS FIRMWARE
    Probe_Selection.value=(Probe_Selection.options[(get_value(table.rows,'name','PCFG_PRB_SEL','value'))])

    tau.value = (np.round((1E6*(2 * math.pi) / ((math.log((((get_value(table.rows,'name', 'PCFG_C11','value')) - 0.5 )/ 16384)))/(-2.71072)*digitizer.ADC_sampling_freq)),1))
    pz.value = np.round( 1E6* (32768 / (digitizer.ADC_sampling_freq*(get_value(table.rows,'name', 'PCFG_CPZ','value'))  )) , 1)
    FLT_CFG_bits=bin((get_value(table.rows,'name', 'PCFG_FLT_CFG','value')))[2:]
    Att1.value=int(str(FLT_CFG_bits[9:13]), 2)
    Att2.value=int(str(FLT_CFG_bits[5:9]), 2)
    BL_IN_Checkbox.value=bool(FLT_CFG_bits[4])
    PZ_Checkbox.value = bool(FLT_CFG_bits[3])
    Filter1_Checkbox.value= bool(FLT_CFG_bits[2])
    Filter2_Checkbox.value= bool(FLT_CFG_bits[1])
    BL_OUT_Checkbox.value= bool(FLT_CFG_bits[0])
    ui.update
    return

def registers_table_update():

    def set_value(rows, key_search, value_search, key_to_modify, new_value):
        """
        Retourne la valeur associée à `key_return` dans la première ligne où
        `key_search` == `value_search`. Retourne None si aucune ligne ne correspond.
        """
        for row in rows:
            if row.get(key_search) == value_search:
                match row['type']:
                    case 'Decimal':
                        row[key_to_modify]= int(new_value)
                    case 'Hexadecimal':
                        row[key_to_modify]= hex(new_value)
        return None

    set_value(table.rows,'name', 'PCFG_THRS','value',(Trigger_Thresh.value))
    #set_value(table.rows,'name', 'PCFG_THRS','value',round(Trigger_Thresh.value/mV_per_ADU))
    set_value(table.rows,'name',  'PCFG_TRG_SOURCE','value',(Trigger_Source.options.index(Trigger_Source.value)))
    set_value(table.rows,'name',  'PCFG_TRG_MODE','value',(Trigger_Mode.options.index(Trigger_Mode.value)))
    set_value(table.rows,'name', 'PCFG_HIST','value',(Trigger_hysteresis.value))
    set_value(table.rows,'name', 'PCFG_TRG_HOLD','value',round((Trigger_hold_off.value)/(1E6/digitizer.ADC_sampling_freq)))
    set_value(table.rows,'name', 'PCFG_TRG_W','value',round(Trigger_output_width.value/(1E6/digitizer.ADC_sampling_freq)))
    set_value(table.rows,'name', 'PCFG_GATE_W','value',round(Gate_window.value/(1E6/digitizer.ADC_sampling_freq)))
    set_value(table.rows,'name', 'PCFG_POL','value',(polarity.options.index(polarity.value)))
    #set_value(table.rows,'name', 'PCFG_RAWOFS','value',(offset.value))
    set_value(table.rows,'name', 'PCFG_RAWOFS','value',round((offset.value)/digitizer.mV_per_ADU))
    set_value(table.rows,'name', 'PCFG_PRB_SEL','value',(Probe_Selection.options.index(Probe_Selection.value)))
    
    TS = (2 * math.pi) / (tau.value*1E-6*digitizer.ADC_sampling_freq)
    pzcoeff2_d =int (32768 / (digitizer.ADC_sampling_freq*pz.value*1E-6))
    coeff11_d = int(math.exp(-2.71072*TS) * 16384 + 0.5)
    coeff12_d = int(2 * math.exp(-1.35536*TS)*math.cos(0.327948*TS) * 16384 + 0.5)
    coeff21_d = int(math.exp(-2.36216*TS) * 16384 + 0.5)
    coeff22_d = int(2 * math.exp(-1.18108*TS)* math.cos(1.06037*TS) * 16384 + 0.5)


    pzcoeff2_h= hex(int(pzcoeff2_d)) 
    coeff11_h=hex(int(coeff11_d))
    coeff12_h=hex(int(coeff12_d))
    coeff21_h=hex(int(coeff21_d))
    coeff22_h=hex(int(coeff22_d))

    set_value(table.rows,'name', 'PCFG_CPZ','value',pzcoeff2_d)
    set_value(table.rows,'name', 'PCFG_C11','value',coeff11_d)
    set_value(table.rows,'name', 'PCFG_C12','value',coeff12_d)
    set_value(table.rows,'name', 'PCFG_C21','value',coeff21_d)
    set_value(table.rows,'name', 'PCFG_C22','value',coeff22_d)

    # set_value(table.rows,'name', 'PCFG_CPZ','value',pzcoeff2_h)
    # set_value(table.rows,'name', 'PCFG_C11','value',coeff11_h)
    # set_value(table.rows,'name', 'PCFG_C12','value',coeff12_h)
    # set_value(table.rows,'name', 'PCFG_C21','value',coeff21_h)
    # set_value(table.rows,'name', 'PCFG_C22','value',coeff22_h)
    
    # FLT_CFG register value for the configuration of the filter stage
    Att1_bits=np.binary_repr(int(Att1.value), width=4)
    Att2_bits=np.binary_repr(int(Att2.value), width=4)
    BL_in_bit=np.binary_repr(BL_IN_Checkbox.value, width=1)
    PZ_bit=np.binary_repr(PZ_Checkbox.value, width=1)
    G1_bit=np.binary_repr(Filter1_Checkbox.value, width=1)
    G2_bit=np.binary_repr(Filter2_Checkbox.value, width=1)
    BL_out_bit=np.binary_repr(BL_OUT_Checkbox.value, width=1)

    FLT_CFG_bits='000'+BL_out_bit+G1_bit+G2_bit+PZ_bit+BL_in_bit+Att2_bits+Att1_bits

    def bin_to_hex(binary_str):
        return hex(int(binary_str, 2))[2:] 

    set_value(table.rows,'name', 'PCFG_FLT_CFG','value',int(bin_to_hex(FLT_CFG_bits),16))

    ui.update(table)
    ui.update()

async def close_program():

    print("Stopping acquisition")
    live_mode_switch.value=False
    await startToggleButton.toggle()

    print("Disconnecting SciSDK")
    digitizer.disconnect()

    # Give time for SSH command to be transmitted
    print("Waiting before closing GUI")
    await asyncio.sleep(1)

    print("Closing browser window")
    await ui.run_javascript('window.close()')

    await asyncio.sleep(0.5)

    print("Stopping NiceGUI server")
    app.shutdown()

######################################################################
# Custom GUI buttons
######################################################################

class noise_ToggleButton(ui.button):
    def __init__(self, *args, **kwargs) -> None:
        
        self._state = False
        super().__init__(*args, **kwargs)
        self.on('click', self.toggle)

    async def toggle(self) -> None:
        """Toggle the button state."""
        if not self._state:
            self._state = not self._state
            self.update()
            await asyncio.sleep(1)
            await start_noise_measurement()
            if self._state:
                self._state = not self._state
                self.set_visibility(1)
            self.update()
            self.enable
        else :
            #self.set_visibility(0)
            self._state = not self._state
            self.update()
            self.disable
            await asyncio.sleep(1)
            #self.set_visibility(1)
            self.enable
            #self.update()  
        #time.sleep(2)

    def update(self) -> None:
        self.props(f'color={"red" if self._state else "blue"}')
        self.props(f'icon={"block" if self._state else "ads_click"}')
        if self._state : 
            self.text="Stop"
        else :
            self.text="Start"
        super().update()

class acqToggleButton(ui.button):
    def __init__(self, *args, **kwargs) -> None:
        
        self._state = False
        super().__init__(*args, **kwargs)
        self.on('click', self.toggle)

    async def toggle(self) -> None:
        """Toggle the button state."""
        if not self._state:
            self._state = not self._state
            self.update()
            #await asyncio.sleep(1)
            await start_acquisition()
            if self._state:
                self._state = not self._state
                self.set_visibility(1)
            self.update()
            self.enable
        else :
            #self.set_visibility(0)
            self._state = not self._state
            self.update()
            self.disable
            await asyncio.sleep(1)
            #self.set_visibility(1)
            self.enable
            #self.update()  
        #time.sleep(2)

    def update(self) -> None:
        self.props(f'color={"red" if self._state else "blue"}')
        self.props(f'icon={"block" if self._state else "ads_click"}')
        if self._state : 
            self.text="Stop"
        else :
            self.text="Start"
        super().update()

####################################################################################
# GUI layout
####################################################################################

# GUI browser tab title
if config.debug_mode:
    ui.page_title('QDIV-gui - Test Mode')
else:
    ui.page_title('QDIV-gui')

# GUI top menu
with ui.row().classes('w-full no-wrap'):
    with ui.row().classes('w-full'):
        
        # GUI title
        with html.section().style('font-size: 250%').classes('no-wrap'):
            html.strong('QDIV-gui ').classes('w-1/8') \
                .classes('cursor-pointer') \
                .tooltip('Python Graphical User Interface (version ' + version +  ') for Position-Sensitive Neutron Detectors read-out by the resistive-charge-division method implemented in the QDIV firmware initially developped at the ILL and running in CAEN digitizers.')
            html.hr()

        with ui.row().classes('gap-6 no wrap '):
            logo_1=ui.image('ILL_logo.png').classes('w-14 h-auto  object-contain top-3')
            logo_2=ui.image('ESS_logo.png').classes('w-20 h-auto  object-contain top-3')
            logo_3=ui.image('ISIS_logo.png').classes('w-20 h-auto  object-contain top-3')
            logo_4=ui.image('FRMII_logo.png').classes('w-20 h-auto  object-contain top-3')

        

# GUI tab selection
    with ui.row().classes('justify-end'):
        with ui.tabs() as tabs:
            ui.tab('a', label='Acquisition', icon='image')
            #ui.tab('s', label='Spectrum', icon='filter')
            ui.tab('o', label='Oscilloscope', icon='filter')
            ui.tab('n', label='Elec. noise', icon='waves').classes("justify-end")
            ui.tab('g', label='Gain uniformity', icon='equalizer').classes("justify-end")
            ui.tab('d', label='Settings', icon='settings')
            ui.tab('X', label='Close', icon='close').classes("justify-end")

# GUI tab panels
with ui.tab_panels(tabs,value='a').classes('w-full'):
    # GUI tab "Acquisition"
    with ui.tab_panel('a'):
        with ui.row().classes('w-full no-wrap'):
            with ui.column().classes('w-1/5 no-wrap'):
                # Acquisition parameters
                live_mode_switch = ui.switch('LIVE MODE',on_change=live_mode)
                acq_Trigger_Thresh=ui.number(label='Filter threshold [0-'+str(2**digitizer.filter_output_bit_depth)+']', value=0, format='%d').classes('w-full no-wrap')
                startToggleButton=acqToggleButton(text='Start').classes('w-full')
                expT=ui.number(label='Acquisition time (s)', value=1.0, format='%.2f',).classes('w-full no-wrap')
                repetitions=ui.number(label='Number of acquisitions', value=1, format='%d').classes('w-full no-wrap')
                
                # File saving parameters
                def plots_saving_click(): 
                    match plots_saving_checkbox.value:
                        case 0:
                            path.visible=False
                            subpath.visible=False
                            list_saving_checkbox.visible=False
                        case 1: 
                            path.visible=True
                            subpath.visible=True
                            list_saving_checkbox.visible=True
                    ui.update
                    return
                
                plots_saving_checkbox = ui.checkbox('Save plots',value=False,on_change=plots_saving_click)
                list_saving_checkbox = ui.checkbox('Save LIST data',value=False).bind_visibility_from(plots_saving_checkbox,'value').classes('w-full')
                subpath=ui.input(label='Subpath',value=config.det.name + '_' +str(config.det.SN).zfill(2)  +'/' + str(date.today())+'/run1',
                validation={'Input too long': lambda value: len(value) < 50}).bind_visibility_from(plots_saving_checkbox,'value').classes('w-full')
                path=ui.input(label='Path',value=config.path,
                validation={'Input too long': lambda value: len(value) < 50}).bind_visibility_from(plots_saving_checkbox,'value').classes('w-full')
        
            with ui.column().classes('w-full justify-center no-wrap items-center'):
                
                acq_fig_selec_toggle = ui.toggle(['Image','Pulse Height Spectra','Single tube','Pulse height vs position'], value='Image',on_change=update_acquisition_figures).classes('justify-center no-wrap')
                with ui.row().classes('w-full justify-end no-wrap items-center'): 
                    tube_number=ui.number(label='Tube number', min=1,max=config.det.nb_of_tubes,value=1, format='%d',on_change=update_acquisition_figures).classes('w-1/8 justify-center no-wrap').bind_visibility_from(acq_fig_selec_toggle, 'value', lambda v: 'Single tube' in v)
                    projection_checkbox = ui.checkbox('Projection',value=False,on_change=update_acquisition_figures).bind_visibility_from(acq_fig_selec_toggle, 'value', lambda v: 'Pulse Height Spectra' in v)
                    log_checkbox=ui.checkbox('log',value=False,on_change=update_acquisition_figures)
                fig=px.imshow(np.ones([config.det.nb_of_tubes,config.nb_of_pixels_per_tube]),aspect=config.det.aspect_ratio,labels=dict(x="X channel", y="Y channel", color="Counts"))
                
                with ui.column().classes('w-full justify-center no-wrap items-center'):

                    with ui.row().classes('w-full justify-center no-wrap items-center'):       
                        ACQ_FIG1=ui.plotly(fig).classes('w-full wrap justify-center').bind_visibility_from(acq_fig_selec_toggle, 'value', lambda v: 'Image' or 'Pulse Height Spectra' in v)
                        #ACQ_FIG2=ui.plotly(fig).classes('w-full justify-center wrap').bind_visibility_from(acq_fig_selec_toggle, 'value', lambda v: 'Image'  in v)
                        
                        ACQ_PROJ=ui.matplotlib(figsize=(3, 2)).bind_visibility_from(acq_fig_selec_toggle, 'value', lambda v: 'Pulse Height Spectra' in v).classes('w-1/2 justify-center no-wrap')                   
                        ACQ_FIG2 = ui.line_plot(n=1, limit=config.nb_of_pixels_per_tube,update_every=config.GUI_refresh_period)\
                                .with_legend(['Ch A', 'Ch B'], loc='upper right', ncol=1).bind_visibility_from(acq_fig_selec_toggle, 'value', lambda v: 'Single tube' in v).classes('w-1/2 no-wrap items-center justify-center')
                        ACQ_FIG2.fig.suptitle('Counts')
                        ACQ_FIG2.fig.supxlabel("Tube length (a.u.)")
                        ACQ_FIG2.fig.supylabel('Counts')
                        
                        ACQ_FIG3 = ui.line_plot(n=3, limit=config.nb_of_bins_in_spectrum,update_every=config.GUI_refresh_period)\
                                .with_legend(['Ch A + Ch B', 'Ch A','Ch B'], loc='upper right', ncol=1).bind_visibility_from(acq_fig_selec_toggle, 'value', lambda v: 'Single tube' in v).classes('w-1/2 no-wrap items-center justify-center')
                        ACQ_FIG3.fig.suptitle('Pulse Height Spectra')
                        ACQ_FIG3.fig.supxlabel("Bins")
                        ACQ_FIG3.fig.supylabel("Counts")

    # GUI tab "Oscilloscope"
    with ui.tab_panel('o'):
        with ui.row().classes('w-full no-wrap'):
            with ui.column().classes('w-1/8 items-center '):

                html.b('Input stage')
                polarity = ui.select(label='Polarity',options=["Positive","Negative"], value="Negative").classes('w-full')
                #offset=ui.number(label='Offset [0-65535]', value=20000, format='%d').classes('w-full')
                offset=ui.number(label='Offset [mV]', value=0, format='%.2f').classes('w-full')
                html.b('Digital Filter')
                Trigger_Thresh=ui.number(label='Filter threshold [0-'+str(2**digitizer.filter_output_bit_depth)+']', value=0, format='%d').classes('w-full')
                #Trigger_Thresh=ui.number(label='Filter threshold [mV]', value=0, format='%d').classes('w-full')
                Trigger_Thresh.bind_value(acq_Trigger_Thresh)

                Trigger_Source = ui.select(label='Source',options=["Filter output","Filter input"], value="Filter output").classes('w-full')
                Trigger_Mode = ui.select(label='Mode',options=["A+B","No internal trigger","A or B","A & B"], value="A+B").classes('w-full')
                Trigger_hysteresis = ui.number(label='Hysteresis [0-65535]', value=1000, format='%d').classes('w-full') 
                Trigger_hold_off = ui.number(label='Hold-off time (us)', value=1, format='%.2f').classes('w-full')
                Trigger_output_width = ui.number(label='Trigger output width (us)', value=0.1, format='%.2f').classes('w-full')
                Gate_window = ui.number(label='Gate window (us)', value=1, format='%.2f').classes('w-full')


            with ui.column().classes('w-1/8 items-center '):
                #html.b('Digital filter')
                registers_table_update_button=ui.button(text='Apply filter settings',on_click=lambda:(registers_table_update(),digitizer.set_registers(table))).classes('w-full')
                tau=ui.number(label='Shaping time (us)', value=1.0, format='%.1f').classes('w-full')
                pz=ui.number(label='Pole Zero (us)', value=1.0, format='%.1f').classes('w-full')
                Att1=ui.number(label='Filter 1 att. [0-7]', value=1, format='%d').classes('w-full')
                Att2=ui.number(label='Filter 2 att. [0-7]', value=1, format='%d').classes('w-full')
                BL_IN_Checkbox = ui.checkbox('BL IN',value=True,on_change=plots_saving_click).classes('w-full')
                PZ_Checkbox = ui.checkbox('PZ',value=True,on_change=plots_saving_click).classes('w-full')
                Filter1_Checkbox = ui.checkbox('Filter 1',value=True,on_change=plots_saving_click).classes('w-full')
                Filter2_Checkbox = ui.checkbox('Filter 2',value=True,on_change=plots_saving_click).classes('w-full')
                BL_OUT_Checkbox = ui.checkbox('BL OUT',value=True,on_change=plots_saving_click).classes('w-full')
            
            with ui.column().classes('w-1/8 items-center '):
                html.b('Oscilloscope')
                
                 
                line_updates = ui.timer(0.1, update_scope_traces, active=False)
                scope_switch = ui.switch('Run / Stop',on_change=lambda:digitizer.reset_pulse_height_spectrum(tube_number.value)).bind_value_to(line_updates,'active')
                ui.label('Scope trigger threshold [mV]')
                #scope_threshold_slider = ui.slider(min=0, max=2**16, value=14000)
                scope_threshold_slider = ui.slider(min=-1000, max=1000, value=-500)
                ui.label().bind_text_from(scope_threshold_slider, 'value')
                scope_trigger_mode=ui.select(label='Scope trigger mode',options=["disabled", "self", "analog", "ext", "digital" ], value="analog").classes('w-full')
                #scope_trigger_edge=ui.select(label='Scoppe trigger edge',options=["Rising", "Falling" ], value="Rising").classes('w-full')
                scope_trig_ch_button =  ui.select(label='Scope trigger channel',options=["Channel A","Channel B"], value="Channel A").classes('w-full')
                scope_tube_number=ui.number(label='Tube number', min=1,max=config.det.nb_of_tubes,value=1, format='%d').classes('w-full').bind_value(tube_number)
                an_mon_sel_options=["signal with corrected offset and polarity",
                "signal with corrected offset only",
                "main filter output",
                "sum of both side main filter output",
                "peak detector output",
                "peak detector output on sum signal",
                "input baseline correction algorithm output",
                "pz algorithm output",
                "first filter algorithm output (divided by 2)",
                "second filter algorithm output (divided by 2)",
                "output baseline algorithm output (divided by 2)"]
                Probe_Selection = ui.select(label='Probe selection',options=an_mon_sel_options, value=an_mon_sel_options[2]).classes('w-full')
                decimation_factor =ui.number(label='Decimation fator', value=0, format='%d').classes('w-full')
                
                # File saving parameters
                def traces_autosave_click(): 
                    match traces_autosave_checkbox.value:
                        case 0:
                            traces_path.visible=False
                            traces_subpath.visible=False
                            number_of_traces_input.visible=False
                        case 1: 
                            traces_path.visible=True
                            traces_subpath.visible=True
                            number_of_traces_input.visible=True
                    ui.update
                    return
                
                traces_autosave_checkbox = ui.checkbox('Save traces',value=False,on_change=traces_autosave_click)
                number_of_traces_input = ui.number(label='Number of traces', value=10, format='%d').classes('w-full')
                traces_subpath=ui.input(label='Subpath',value=config.det.name + '_' +str(config.det.SN).zfill(2)  +'/' + str(date.today())+'/traces/',
                validation={'Input too long': lambda value: len(value) < 50})
                traces_path=ui.input(label='Path',value=config.path,
                validation={'Input too long': lambda value: len(value) < 50})
                number_of_traces_input.visible=False
                traces_path.visible=False
                traces_subpath.visible=False
            
            #with ui.column().classes('w-3/8 justify-center no-wrap'): 
            with ui.column().classes(' no-wrap items-center justify-center'): 
                stats_or_spectra_toggle = ui.toggle(['Filter input/output','Digital signals','Statistics', 'A+B pulse height spectrum'], value='Filter input/output').classes('justify-center no-wrap')
                
                with ui.row().classes('w-full no-wrap'):

                    with ui.column().classes('w-full no-wrap items-center justify-center'):

                        analog_input_line_plot = ui.line_plot(n=2, limit=digitizer.custom_packet_buffer_size,update_every=config.GUI_refresh_period)\
                            .with_legend(['Ch A', 'Ch B'], loc='upper right', ncol=1)
                        analog_input_line_plot.fig.suptitle("Filter input signals (Tube "+str(tube_number.value)+ ')')
                        analog_input_line_plot.fig.supxlabel("Micro-seconds")
                        analog_input_line_plot.fig.supylabel("mV")

                        with ui.row().classes('no-wrap items-center justify-center').bind_visibility_from(stats_or_spectra_toggle, 'value', lambda v: 'Statistics' in v):
                            
                            with ui.column().classes('w-1/2 no-wrap items-center justify-center'): 
                                html.b('Noise [mV rms]')
                                rms_val_ch1=ui.number(label='Channel A input :', value=0.0, format='%.2f').classes('w-full')
                                rms_val_ch2=ui.number(label='Channel B input :', value=0.0, format='%.2f').classes('w-full')
                            
                            with ui.column().classes('w-1/2 no-wrap items-center justify-center'): 
                                html.b('Offset [mV]')
                                mean_val_ch1=ui.number(label='Channel A input :', value=0.0, format='%.2f').classes('w-full')
                                mean_val_ch2=ui.number(label='Channel B input :', value=0.0, format='%.2f').classes('w-full')

                    with ui.column().classes('w-full no-wrap items-center justify-center'):

                        analog_output_line_plot = ui.line_plot(n=2, limit=digitizer.custom_packet_buffer_size,update_every=config.GUI_refresh_period)\
                            .with_legend(['Ch A', 'Ch B'], loc='upper right', ncol=1)
                        analog_output_line_plot.fig.suptitle("Filter output signals (Tube "+str(tube_number.value)+ ')')
                        analog_output_line_plot.fig.supxlabel("Micro-seconds")
                        analog_output_line_plot.fig.supylabel("ADU")

                        with ui.row().classes('no-wrap items-center justify-center').bind_visibility_from(stats_or_spectra_toggle, 'value', lambda v: 'Statistics' in v):
                            
                            with ui.column().classes('w-1/2 no-wrap items-center justify-center'): 

                                html.b('Noise [ADUs rms]')
                                rms_val_ch3=ui.number(label='Channel A output :', value=0.0, format='%.2f').classes('w-full')
                                rms_val_ch4=ui.number(label='Channel B output :', value=0.0, format='%.2f').classes('w-full')
                            
                            with ui.column().classes('w-1/2 no-wrap items-center justify-center'): 

                                html.b('Offset [ADUs]')
                                mean_val_ch3=ui.number(label='Channel A output :', value=0.0, format='%.2f').classes('w-full')
                                mean_val_ch4=ui.number(label='Channel B ouput :', value=0.0, format='%.2f').classes('w-full')
                            
                        digital_line_plot = ui.line_plot(n=4, limit=digitizer.custom_packet_buffer_size, update_every=config.GUI_refresh_period) \
                            .with_legend(['Filter trigger', 'Scope over thresh.','Gate','Energy ready'], loc='upper right', ncol=1).bind_visibility_from(stats_or_spectra_toggle, 'value', lambda v: 'Digital signals' in v)
                        digital_line_plot.fig.suptitle("Digital signals")
                        digital_line_plot.fig.suptitle("Digital signals")
                        digital_line_plot.fig.supxlabel("Micro-seconds")
                        digital_line_plot.fig.supylabel("0/1")

                        reset_spectra_button=ui.button(text='Reset spectra',on_click=lambda:digitizer.reset_pulse_height_spectrum(tube_number.value)).bind_visibility_from(stats_or_spectra_toggle, 'value', lambda v: 'A+B pulse height spectrum' in v)
                        spectra_plot = ui.line_plot(n=1, limit=digitizer.custom_packet_buffer_size,update_every=config.GUI_refresh_period)\
                            .with_legend(['Ch A + Ch B']).classes('w-full justify-center no-wrap').bind_visibility_from(stats_or_spectra_toggle, 'value', lambda v: 'A+B pulse height spectrum' in v)
                        spectra_plot.fig.suptitle(" Pulse height spectrum (Tube "+str(tube_number.value)+ ')')
                        spectra_plot.fig.supxlabel("ADU")
                        spectra_plot.fig.supylabel("Bins")

                
    # GUI menu "Electronic noise"
    with ui.tab_panel('n'):
        with ui.row().classes('w-full no-wrap'):
            with ui.column().classes('w-1/3'):
                    noise_meas_button=noise_ToggleButton(text='Start').classes('w-full')
                    number_of_traces_per_tube_input = ui.number(label='Number of scope traces per tube', value=5, format='%d').classes('w-full').bind_value(number_of_traces_input, 'value') 
                    noise_subpath=ui.input(label='Subpath',value=config.det.name + '_' +str(config.det.SN).zfill(2)  +'/' + str(date.today())+'/elec_noise/',
                    validation={'Input too long': lambda value: len(value) < 50}).classes('w-full')
                    noise_path=ui.input(label='Path',
                    validation={'Input too long': lambda value: len(value) < 50}).bind_value(path,'value').classes('w-full')
            with ui.column().classes('w-full justify-center no-wrap items-center'):
                fig_noise_and_offset_vin=[px.line(np.zeros([config.det.nb_of_tubes,2]),markers=True),px.line(np.zeros([config.det.nb_of_tubes,2]),markers=True)] 
                fig_noise_and_offset_vout=fig_noise_and_offset_vin
                noise_or_offset_toggle = ui.toggle(['RMS noise','Offsets'], value='RMS noise').classes('justify-center no-wrap')
                with ui.row().classes('w-full justify-center no-wrap items-center'):
                    #FIG_NOISE_AND_OFFSET=ui.plotly(fig_noise_and_offset[0]).classes('w-full justify-center no-wrap')
                    FIG_NOISE_AND_OFFSET_VIN=ui.plotly(fig_noise_and_offset_vin[0]).classes('w-1/2 justify-center items-center no-wrap')
                    FIG_NOISE_AND_OFFSET_VOUT=ui.plotly(fig_noise_and_offset_vout[0]).classes('w-1/2 justify-center items-center no-wrap')
                #p=ui.pagination(1,2,direction_links=True,on_change=update_fig_noise_and_offset).classes('w-full justify-center no-wrap')
                FIG_NOISE_AND_OFFSET_VIN.set_visibility(True)
                FIG_NOISE_AND_OFFSET_VOUT.set_visibility(True)

    # # GUI menu "Gain uniformity"
    # #    with ui.tab_panel('g'):
    #     # to be completed

    # GUI menu "Digitizer"
    with ui.tab_panel('d'):


        with ui.row().classes('w-full no-wrap centre-items justify-center'):
            with ui.column().classes('w-1/4 no-wrap justify-left'):
                html.b('Detector')
                # Detector type selection
                def det_type_change():
                    config.det= next((det for det in detectors.list if det.name == det_type.value), None)
                    ui.notify(str(det_type.value) + " selected", timeout=0.2)
                    digitizer.test()
                    return     
                det_type = ui.select([det.name for det in detectors.list],label='Detector',value=config.det.name,on_change=det_type_change).classes('w-full no-wrap justify-center')
                
                # Detector serial number selection
                def SN_change():
                    ui.notify("Serial Number modified to " + str(serial_number.value), timeout=0.2)
                    return  
                serial_number = ui.number(label='Serial Number',value=config.det.SN,on_change=SN_change).classes('w-full  no-wrap justify-center')


                html.b('Digitizer')
                shared_values = {"shared_value": "Initial Value"}
                hardware_selec=ui.select(label='Model',options=["R5560","dt1260"], value=config.hardware).classes('w-full  no-wrap justify-center')
                firmware_selec=ui.select(label='Firmware',options=["cspec_rmm 2026.03.11"], value=config.firmware).classes('w-full  no-wrap justify-center')
                regs_def_file_selec=ui.input(label='Registers definition filename', value=config.registers_def_file).classes('w-full  no-wrap justify-center')

            with ui.column().classes('w-1/4  no-wrap justify-center'):
                register_filename=ui.input(label='Registers settings filename', value=config.register_filename,
                validation={'Input too long': lambda value: len(value) < 50}).classes('w-full')

                def openRegisterFile_button_on_click():
                    reg_rows,reg_columns=openRegisterFile(register_filename.value)
                    refresh_table(reg_rows)
                    ui.notification("Register file "+ register_filename.value + " opened !",timeout=1)

                def saveRegisterFile():
                    table.update()
                    reg_items = ['name', 'type', 'address', 'value']
                    with open(register_filename.value, 'w', newline='') as csvfile:
                        writer = csv.writer(csvfile,delimiter=';')
                        writer.writerows([reg_items])
                        [writer.writerows([table.rows[n].values()]) for n in range(np.size(table.rows))]
                    ui.notification("Register file "+ register_filename.value + " saved !",timeout=1)
                
                open_regs_file_button=ui.button(text='open register file',on_click=lambda:openRegisterFile_button_on_click()).classes('w-full')
                save_regs_file_button=ui.button(text='save register file',color='red',on_click=lambda:saveRegisterFile()).classes('w-full')


                #port_selec=ui.select(label='Port',options=["eth","usb"], value="eth").classes('w-full')
                # match hardware_selec.value:
                #     case "R5560":
                #         number_of_ports=round(config.det.nb_of_tubes/16)
                #         IP_addr=[]
                #         for port_number in range(number_of_ports):
                #             IP_addr.append(ui.input(label='IP address ' + str(port_number), value=config.IP[port_number]).classes('w-1/3'))

            with ui.column().classes('w-1/2 no-wrap justify-center'):
                html.b('Registers settings')
                # Function to refresh the table after edits
                def refresh_table(reg_rows):
                    table.rows = reg_rows

                # Function to handle inline edits
                def update_cell(row_id, field, value):
                    for r in reg_rows:
                        if r['id'] == row_id:
                            r[field] = value
                            break
                    ui.notify(f"Updated row {row_id}: {field} = {value}")


                def openRegisterFile(register_filename):
                    reg_dict={};reg_rows=[];reg_items = ['name', 'type', 'address', 'value']
                    reg_columns= [
                        {'name': 'name', 'label': 'Name', 'field': 'name', 'align': 'left'},
                        {'name': 'type', 'label': 'Type', 'field': 'type', 'align': 'right'},
                        {'name': 'address', 'label': 'Address', 'field': 'address', 'align': 'right'},
                        {'name': 'value', 'label': 'Value', 'field': 'value', 'align': 'right'},
                    ]
                    try:
                        with open(register_filename, 'r') as file:
                            next(file)
                            reader = csv.reader(file, delimiter=';')
                            for row in reader:
                                if len(row) == 5:
                                    row.pop(3)  # L'index 3 correspond au quatrième élément 
                                reg_dict={reg_items[index]:value for index,value in enumerate(row)}
                                reg_rows.append(reg_dict)
                            
                    except FileNotFoundError:
                        print("Error: File not found")
                    
                    return reg_rows,reg_columns

                # Create the table
                (reg_rows,reg_columns)=openRegisterFile(config.register_filename)
                table = ui.table(columns=reg_columns, rows=reg_rows, row_key='id', pagination=6).classes('w-full')
                

                # Editable Name column
                table.add_slot('body-cell-name', r'''
                <q-td :props="props">
                    <q-input dense v-model="props.row.name" @update:model-value="val => $parent.$emit('cell-edit', {id: props.row.id, field: 'name', value: val})" />
                </q-td>
                ''')

                # Editable Age column
                table.add_slot('body-cell-age', r'''
                <q-td :props="props">
                    <q-input type="number" dense v-model.number="props.row.age" @update:model-value="val => $parent.$emit('cell-edit', {id: props.row.id, field: 'age', value: val})" />
                </q-td>
                ''')

                # Event listeners for edits and deletes
                table.on('cell-edit', lambda e: update_cell(e.args['id'], e.args['field'], e.args['value']))
            

    # GUI menu "close"
    with ui.tab_panel('X'):
        ui.label('Close the GUI and exit the program')
        ui.button('Close', on_click=close_program).classes('w-1/8')


####################################################################################
# GUI start-up
####################################################################################
digitizer.set_registers(table)
update_GUI_scope_settings()
scope_settings={'tube_number':tube_number.value,'threshold_slider':scope_threshold_slider.value,'trigger_mode':scope_trigger_mode.value,'trig_ch': scope_trig_ch_button.options.index(scope_trig_ch_button.value),'decimation_factor':decimation_factor.value}
ui.run(reload=False)
