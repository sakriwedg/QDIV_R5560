# -*- coding: utf-8 -*-
"""
CAEN digitizers control functions for QDIV-gui

Created on Sat Jul 11 07:44:11 2026

@author: marchalj (ILL)
"""
version="11.07.2026"

from scisdk.scisdk import SciSDK
from scisdk.scisdk_defines import *
import math
import numpy as np
import time
import subprocess


import config # python file with QDIV_gui configuration parameters (file located in QDIV_gui folder) 

# Digitizer class
class digitizer:
    """
    Defines CAEN digitizers parameters and functions (depending on hardware and firmware)

    Attributes
    ----------
    none
        --> digitizer hardware and firmware parameters are defined in config.py

    Methods
    -------
    set_registers(self,registers_table)
    read_scope_traces(self,scope_settings)
    read_pulse_height_spectrum(self,tube_number)
    reset_pulse_height_spectrum(self,tube_number)
    configure_data_readout(self):
    def start(self)
    def stop(self)
    def read_data(self)

    """
    def __init__(self):            
        match config.hardware:
            case "R5560":
                self.ADC_sampling_freq=120E6 #MHz
                self.ADC_input_range=2 #V
                self.ADC_bit_depth=14
                self.filter_output_bit_depth=16
                self.mV_per_ADU=self.ADC_input_range*1000/2**self.ADC_bit_depth
                self.custom_packet_buffer_size=1024
                self.sdk_list=[]
                self.custom_packet_buffer_list=[]
                self.custom_packet = "CP_ETHERNET"
                for i,RJ45_port in enumerate(config.IP):
                    sdk = SciSDK()  # initialize scisdk library
                    # "board0" (registers, oscilloscope, spectra) and "board0list" (event list
                    # read-out) are two device names of the SAME sdk instance : the instance must
                    # be appended only once, sdk_list holds one entry per RJ45 port (it is indexed
                    # by port number in set_registers, read_scope_traces and read_data)
                    res = sdk.AddNewDevice(RJ45_port+':8888', "R5560",config.registers_def_file, "board0")
                    res = sdk.AddNewDevice(RJ45_port+':8888', "R5560",config.registers_def_file, "board0list")
                    self.sdk_list.append(sdk)

                    print(self.sdk_list)
                if res != 0:
                    print(" ! Failed to connect, code:", res)
                    exit(1)
                if config.chatty:
                    library_version = sdk.GetLibraryVersion()
                    print("R5560 library version: " + str(library_version))
                    print("R5560 firmware version: " + str(config.firmware))
#               TODO: Read firmware version from digitizer itself in the GUI                         
            case "dt1260":
                print("Connect to dt1720")

    def set_registers(self,registers_table):
        match config.hardware:
            case "R5560":
                match config.firmware:
                    case "cspec_rmm 2026.03.11":
                        for i,RJ45_port in enumerate(config.IP):
                            for params in registers_table.rows: 
                                match params["type"]:
                                    case "Decimal":
                                        value_int=int(params["value"])
                                    case "Hexadecimal":
                                        value_int=int(params["value"],16)
                                # Extract register short name (there are 3 prefix cases : EMPTY, PCFG_ and FINE_OFS_)
                                if params["name"].startswith("PCFG_"):
                                    param_name = params["name"].replace("PCFG_", "")
                                    path = f"board0:/MMCComponents/PCFG.{param_name}"
                                    if config.chatty:
                                        print(f"Setting {param_name} at path {path} to value {value_int}")
                                    err = self.sdk_list[i].SetParameterInteger(path, value_int)
                                else:
                                    param_name = params["name"]
                                    path = f"board0:/Registers/{param_name}"
                                    if config.chatty:
                                        print(f"Setting {param_name} at path {path} to value {value_int}")
                                    err = self.sdk_list[i].SetRegister(path, int(value_int))    
                                    if config.chatty:
                                        err, value = self.sdk_list[i].GetRegister(path)
                                        print("Register set to {0}\n".format(value))

    def read_scope_traces(self,scope_settings):
        RJ45_port=math.floor((scope_settings['tube_number']-1)/16)
        sdk=self.sdk_list[RJ45_port]
        tube_number_in_port=int((scope_settings['tube_number']-1)-RJ45_port*16)
        oscilloscope_channel="Oscilloscope_"+str(tube_number_in_port)
        res = sdk.SetParameterString("board0:/MMCComponents/"+oscilloscope_channel+".data_processing","decode")
        res = sdk.SetParameterInteger("board0:/MMCComponents/"+oscilloscope_channel+".trigger_level", round(scope_settings['threshold_slider']/self.mV_per_ADU+(2**self.ADC_bit_depth)/2))
        res = sdk.SetParameterString("board0:/MMCComponents/"+oscilloscope_channel+".trigger_mode",scope_settings['trigger_mode'])
        res = sdk.SetParameterInteger("board0:/MMCComponents/"+oscilloscope_channel+".trigger_channel",scope_settings['trig_ch'])
        res = sdk.SetParameterInteger("board0:/MMCComponents/"+oscilloscope_channel+".pretrigger", 500)
        res = sdk.SetParameterInteger("board0:/MMCComponents/"+oscilloscope_channel+".decimator", int(scope_settings['decimation_factor']))
        res = sdk.SetParameterString("board0:/MMCComponents/"+oscilloscope_channel+".acq_mode", "blocking")
        res = sdk.SetParameterInteger("board0:/MMCComponents/"+oscilloscope_channel+".timeout", 500)
        # allocate buffer for oscilloscope
        res, buf_osc = sdk.AllocateBuffer("board0:/MMCComponents/"+oscilloscope_channel)
        res_osc, buf_osc_local = sdk.ReadData("board0:/MMCComponents/"+oscilloscope_channel, buf_osc)
        if res_osc == 0:
            # -------- ANALOG --------
            samples = buf_osc_local.info.samples_analog
            channels = buf_osc_local.info.channels
            analog_time_axis=(1E6*np.arange(samples)/(self.ADC_sampling_freq))*(scope_settings['decimation_factor']+1)
            analog_traces = np.zeros([channels,samples])
            for ch in range(channels):
                start = ch * samples
                end = start + samples
                data = [buf_osc_local.analog[i] & 0xFFFF for i in range(start, end)]
                analog_traces[ch]=data
            # -------- DIGITAL --------
            samples = buf_osc_local.info.samples_digital
            tracks = buf_osc_local.info.tracks_digital_per_channel
            DIGITAL_CHANNEL=2
            #digital_time_axis=np.arange(samples)/(ADC_SAMPLING_FREQ)
            digital_traces = np.zeros([tracks,samples])
            for d in range(tracks):
                digital_wave = []
                for i in range(samples):
                    index = DIGITAL_CHANNEL * tracks * samples + d * samples + i
                    raw = buf_osc_local.digital[index] & 0xFF
                    bit = raw & 1
                    digital_wave.append(0.5 * bit + (tracks - d - 1))
                digital_traces[d, :] = digital_wave
        else:
            analog_time_axis=[]
            analog_traces=[]
            digital_traces=[]
        return(res_osc,analog_time_axis,analog_traces,digital_traces)

    def read_pulse_height_spectrum(self,tube_number):
        RJ45_port=math.floor((tube_number-1)/16)
        sdk=self.sdk_list[RJ45_port]
        tube_number_in_port=int((tube_number-1)-RJ45_port*16)
        spectrum_channel="Spectrum_"+str(tube_number_in_port)
        # Acquire spectra
        sdk.SetRegister("board0:/Registers/noisepower", 100)
        sdk.SetRegister("board0:/Registers/centroid", 1000)
        # set board parameters
        sdk.SetParameterString("board0:/MMCComponents/" +  spectrum_channel  + ".rebin", "0")
        sdk.SetParameterString("board0:/MMCComponents/"  + spectrum_channel   + ".limitmode", "freerun")
        sdk.SetParameterString("board0:/MMCComponents/" +  spectrum_channel +  ".limit", "100")
        sdk.SetParameterString("board0:/MMCComponents/" +  spectrum_channel  + ".min", "0")
        sdk.SetParameterString("board0:/MMCComponents/" +  spectrum_channel  + ".max", "2**16")
        # execute command start
        sdk.ExecuteCommand("board0:/MMCComponents/" + spectrum_channel +".start", "")
        # allocate buffer
        res, buf = sdk.AllocateBuffer("board0:/MMCComponents/"+spectrum_channel)
        res, buf = sdk.ReadData("board0:/MMCComponents/"+spectrum_channel, buf)# read data from board
        if res == 0:
            bins = []
            phs = []
            for index in range(buf.info.valid_bins):
                bins.append(index)
                phs.append(buf.data[index])
            print(len(phs))
        return(bins,np.array(phs))

    def reset_pulse_height_spectrum(self,tube_number):
        RJ45_port=math.floor((tube_number-1)/16)
        sdk=self.sdk_list[RJ45_port]
        tube_number_in_port=int((tube_number-1)-RJ45_port*16)
        spectrum_channel="Spectrum_"+str(tube_number_in_port)
        sdk.ExecuteCommand("board0:/MMCComponents/" +   spectrum_channel  +".reset", "")
        sdk.ExecuteCommand("board0:/MMCComponents/" + spectrum_channel  +".start", "")
      
    def configure_data_read_out(self):
        # called again at every acquisition start : without this the buffers allocated by the
        # previous acquisitions pile up in the list (and are leaked)
        self.custom_packet_buffer_list=[]
        for i,sdk in enumerate(self.sdk_list):
            sdk.SetParameterString(f"board0list:/MMCComponents/{self.custom_packet}.thread","false")
            sdk.SetParameterString(f"board0list:/MMCComponents/{self.custom_packet}.acq_mode","non-blocking")
            if config.chatty:
                res, v = sdk.GetParameterString(f"board0list:/MMCComponents/{self.custom_packet}.acq_mode")
                print("acq_mode =", v)
            # Allocate CP buffer
            if config.chatty:
                print("\n --- Allocating CP buffer:")
            res, buf_cus = sdk.AllocateBuffer(f"board0list:/MMCComponents/{self.custom_packet}", self.custom_packet_buffer_size*10)
            self.custom_packet_buffer_list.append(buf_cus)
            if config.chatty:
                if res != 0:
                    print("Failed to allocate CP buffer")
                    exit(1)
                else:
                    print(" - CP buffer allocated successfully")

    def start_data_read_out(self):
        for i,sdk in enumerate(self.sdk_list):
            res = sdk.ExecuteCommand(f"board0list:/MMCComponents/{self.custom_packet}.start", "")
        return time.time()

    def stop_data_read_out(self):
        for i,sdk in enumerate(self.sdk_list):
            res = sdk.ExecuteCommand(f"board0list:/MMCComponents/{self.custom_packet}.stop", "")
        return time.time()

    def read_data(self):
        def swap32(x):
            return ((x & 0xFF) << 24) | \
                ((x & 0xFF00) << 8) | \
                ((x & 0xFF0000) >> 8) | \
                ((x >> 24) & 0xFF)
        #timestamp_list=[]
        channel_list=[]
        energy_A_list=[]
        energy_B_list=[]
        read_time_list=[]
        for i,sdk in enumerate(self.sdk_list):
            res_cus, buf = sdk.ReadData(f"board0list:/MMCComponents/{self.custom_packet}",self.custom_packet_buffer_list[i])
            read_time_list.append(time.time())
            valid = int(buf.info.valid_data)
            if config.chatty:
                # never print unconditionally here : this runs in the read-out thread, several
                # hundred times per second, and console output would throttle the read-out
                print('######## Buffer size : '+ str(valid))
            for evt in range(valid):
                #w0 = swap32(buf.data[evt].row[0])
                #w1 = swap32(buf.data[evt].row[1])
                w2 = swap32(buf.data[evt].row[2])
                w3 = swap32(buf.data[evt].row[3])
                #timestamp = (w1 << 32) | w0
                channel   = (w2 >> 8) & 0xFF
                energy_A  = w3 & 0xFFFF
                energy_B  = (w3 >> 16) & 0xFFFF
                # shift channel number to the right tube number (for multiple RJ45 ports)
                channel = channel + 16*i
                #print(f"Event {evt}: Channel {channel}, Energy A {energy_A}, Energy B {energy_B}")
                #timestamp_list.append(timestamp)
                channel_list.append(channel)
                energy_A_list.append(energy_A)
                energy_B_list.append(energy_B)
                
        
        return(channel_list,energy_A_list,energy_B_list,np.mean(read_time_list))


    def disconnect(self):

        def is_reachable(ip):
            try:
                result = subprocess.run(
                    ["ping", "-n", "1", "-w", "1000", ip],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )
                return result.returncode == 0
            except Exception as e:
                print(f"Ping error for {ip}: {e}")
                return False

        def reboot_all_r5560():
            for ip in config.IP:
                print(f"Checking {ip} ...")

                if is_reachable(ip):
                    print(f"{ip} reachable, rebooting")
                    reboot_r5560(ip)
                else:
                    print(f"{ip} not reachable, skipping reboot")


        def reboot_r5560(ip_address="10.128.0.50"):
            """
            Reboot R5560 through SSH.
            The SSH connection will normally close because the device reboots.
            """

            print("\n===================================")
            print(" Reboot procedure started")
            print("===================================")

            ssh_command = [
                "ssh",
                "-o", "HostKeyAlgorithms=+ssh-rsa",
                f"root@{ip_address}",
                "/bin/busybox reboot"
            ]

            print(f"Executing SSH command: {' '.join(ssh_command)}")

            try:
                print("\nSending reboot command...")

                result = subprocess.run(
                    ssh_command,
                    capture_output=True,
                    text=True,
                    timeout=10
                )


                if result.stdout:
                    print("STDOUT:")
                    print(result.stdout)

                if result.stderr:
                    print("STDERR:")
                    print(result.stderr)

                # Normally SSH exits with an error because reboot closes the socket
                if "Connection to" in result.stderr and "closed" in result.stderr:
                    print("R5560 reboot initiated successfully")

                elif result.returncode == 0:
                    print("R5560 reboot command completed")

                else:
                    print("WARNING: reboot command may have failed")


            except subprocess.TimeoutExpired:
                # This is actually a good sign:
                # the device may have rebooted while SSH was waiting
                print("SSH timeout")
                print("Assuming R5560 is rebooting")

            except Exception as e:
                print("Unexpected error during reboot:")
                print(e)

            print("===================================")
            print(" Reboot procedure finished")
            print("===================================\n")

                    

        print("Disconnecting SciSDK")
        for i,sdk in enumerate(self.sdk_list):
            err = sdk.DetachDevice("board0")
            print("DetachDevice returned:", err)

            err = sdk.FreeLib()
            print("FreeLib returned:", err)

        print("Rebooting R5560 boards")
        reboot_all_r5560()


