#!/usr/bin/env python
"""
Sigproc utilities

See here for format description: http://sigproc.sourceforge.net/sigproc.pdf
Copyright (C) Keith Bannister 2011
"""

import struct
import os
import numpy as np
import warnings

global _verbose
_verbose = False
def _debug(msg):
    global _verbose
    if _verbose:
        print(msg)
        
HEADER_LENGTH=1024*10
INT_PARAMS=('telescope_id','machine_id','data_type','barycentric','pulsarcentric','nbits','nsamples', 'nchans', 'nifs')
DOUBLE_PARAMS = ('az_start','za_start','src_raj','src_dej','tstart','tsamp','fch1','foff','fchannel','refdm','period')
STRING_PARAMS = ('rawdatafile','source_name')
FREQ_PARAMS = ('FREQUENCY_START','FREQUENCY_END')
DATA_TYPES = ("UNKNOWN", "filterbank", "time series")
TELESCOPE_IDS = ("fake data", "Arecibo", "Ooty")
MACHINE_IDS = ("FAKE", "PSPM", "WAPP", "OOTY")

INT_FORMAT = 'i'
STRING_FORMAT = 's'
DOUBLE_FORMAT = 'd'

def sigproc_sex2deg(x):
    '''
    Computes a decimal number in degrees from a rediculous sigproc number in the form
    ddmmss.s

    >>> print '%0.5f' % sigproc_sex2deg(123456.789)
    12.58244
    >>> print '%0.5f' % sigproc_sex2deg(-123456.789)
    -12.58244
    >>> sigproc_sex2deg(0)
    0.0
    '''
    sign = 1.
    if x < 0:
        sign = -1.

    x = abs(x)
    
    dd = float(int(x/10000))
    mm = float(int((x - dd*10000)/100))
    ss = x - dd*10000 - mm*100
    dout = sign*(dd + mm/60.0 + ss/3600.0)
    
    return dout

def unpack(hdr, param, struct_format):
    idx = hdr.find(param)
    if idx < 0:
        #warnings.warn('Could not find parameter {}'.format(param))
        return None

    idx += len(param) # idx is thte start of the string
    size = struct.calcsize(struct_format)
    bits = hdr[idx:idx+size]
    value = struct.unpack(struct_format, bits)[0]
    
    #print 'param', param, 'size', size, 'part', bits, 'value', value
    return value

def unpack_str(hdr, param):
    idx = hdr.find(param)
    if idx < 0:
        #warnings.warn('Could not find parameter {}'.format(param))
        return None

    idx += len(param) # idx is thte start of the string
    
    # The size of the string is given in the next field, then the value of the string
    count = struct.unpack('i', hdr[idx:idx+4])[0]
    end_idx = idx+4+count
    value = hdr[idx+4:end_idx]
    
    #print 'param', param, 'start', idx, 'end', end_idx, 'value', value
    return value

def write_str(f, s):
    n = struct.pack('i', len(s))
    f.write(n)
    f.write(s)

def write(f, v, struct_format):
    try:
        f.write(struct.pack(struct_format, v))
    except:
        print('Could not write value %s to %s with format %s' %( v, f, struct_format))
        raise 

class SigprocFile(object):
    def __init__(self, filename, mode='r', header=None):
        self.filename = filename
        self.fin = open(self.filename, mode)
        if header is not None:
            self._write_header(header)
            self.header = header
        else:
            self._read_header()

        if 'src_raj' in self.header:
            self.src_raj_deg = sigproc_sex2deg(self.header['src_raj'])*15.0
            
        if 'src_dej' in self.header:
            self.src_dej_deg = sigproc_sex2deg(self.header['src_dej'])

    def _write_header(self, header):
        f = self.fin
        f.seek(0)
        write_str(f, 'HEADER_START')

        for k, v in header.items():
            if v is None:
                continue
            if k in STRING_PARAMS:
                write_str(f, k)
                write_str(f, v)
            elif k in INT_PARAMS:
                write_str(f, k)
                write(f, v, INT_FORMAT)
            elif k in DOUBLE_PARAMS:
                write_str(f, k)
                write(f, v, DOUBLE_FORMAT)
            else:
                print('Cannot write header', k)

        write_str(f, 'HEADER_END')
        self.data_start_idx = f.tell()
        
    def _read_header(self):
        fin = self.fin
        fin.seek(0)
        hdr = fin.read(HEADER_LENGTH)
        start_idx = hdr.find("HEADER_START")
        end_idx = hdr.find("HEADER_END")
        hdr = hdr[start_idx:end_idx]
        self.data_start_idx = end_idx + len("HEADER_END")
        self.seek_data()
        header = {}
        self.header = header
        self.hdr = hdr
        
                    
        for p in STRING_PARAMS:
            header[p] = unpack_str(hdr, p)
            
        for p in INT_PARAMS:
            header[p] = unpack(hdr, p, INT_FORMAT)
            
        for p in DOUBLE_PARAMS:
            header[p] = unpack(hdr, p, DOUBLE_FORMAT)

        for k,v in self.header.items():
            setattr(self, k, v)
            

        self.file_size_bytes = os.path.getsize(self.filename)
        self.header_size_bytes = self.data_start_idx
        self.data_size_bytes = self.file_size_bytes - self.header_size_bytes
        assert self.nifs > 0, 'Invalid nifs {}'.format(self.nifs)
        assert self.nchans > 0, 'Invalid nchans {}'.format(self.nchans)
        assert self.nbits > 0, 'Invalid nbits {}'.format(self.nbits)
        self.bytes_per_element = self.nifs * self.nchans * self.nbits/8 
        self.file_size_elements = self.data_size_bytes / self.bytes_per_element
        
        if self.nsamples is None:
            self.nsamples = self.file_size_elements
            
        self.observation_duration = self.nsamples * self.tsamp

        
    def seek_data(self, offset_bytes=0):
        self.fin.seek(self.data_start_idx + offset_bytes)

    def seek_sample(self, sampnum):
        self.fin.seek(self.data_start_idx + sampnum*self.nifs*self.nchans*self.nbits/8)
        
    def get_num_elements(self):
        nelements = self.nifs * self.nchans * self.nsamples
        return nelements
    
    def arr_index(self, sampindex, chanindex=0, ifindex=0):
        if chanindex < 0 or chanindex >= self.nchans:
            raise ValueError('Invalid channel index')
        
        if ifindex < 0 or ifindex >= self.nifs:
            raise ValueError('Invalid IF index')
        
        if sampindex < 0 or sampindex >= self.nsamples:
            raise ValueError('Invalid sample index')
        
        arridx = sampindex * self.nifs * self.nchans + ifindex * self.nchans + chanindex
        
        return arridx
    
    def sky_freq(self, chanidx):
        """REturns the sky frequency in  MHz"""
        if chanidx < 0 or chanidx >= self.nchans:
            raise ValueError('Invalid channel index')
        
        freq = self.fch1 + chanidx * self.foff
        
        return freq
    
    def get_data_type(self):
        if self.data_type < len(DATA_TYPES):
            return DATA_TYPES[self.data_type]
        else:
            return None
        
    def get_data(self, time_slice, chanindex=0, ifindex=0):
        if self.nifs != 1:
            raise NotImplementedError("Can't handle Nif > 1")
            
        if self.nbits == 8:
            dtype = np.uint8
            samps_per_element = 1
        elif self.nbits == 32:
            dtype = np.float32
            samps_per_element = 1
        elif self.nbits == 2:
            dtype = np.uint8
            samps_per_element = 4
        else:
            raise NotImplementedError("Can't handle nbits: %d" % self.nbits)
        
        if time_slice.step is not None:
            raise NotImplementedError("Can only handle contiguous slices")
        
        if time_slice.start is None:
            time_start = 0
        else:
            time_start = time_slice.start
            
        if time_slice.stop is None:
            time_end = self.nsamples
        else:
            time_end = time_slice.stop
            
        num_samples = (time_end - time_start)
        num_elements = num_samples*self.nchans
        num_dtypes = num_elements / samps_per_element 
        
        if num_samples < 0:
            raise ValueError("cant do negative number of samples")
        
        byte_start = self.arr_index(time_start, chanindex, ifindex) * self.nbits / 8
        self.seek_data(byte_start)
        
        """ TODO: If nbits < 8 will need to read an extra byte here"""
        num_bytes = num_samples * self.nbits / 8

        if num_bytes + self.data_start_idx > self.file_size_bytes:
            raise ValueError('Requested data off the end of the file. File size: %d. num_bytes: %d. data_start_idx: %d' % (self.file_size_bytes, num_bytes, self.data_start_idx))

        data = np.fromfile(self.fin, dtype=dtype, count=num_dtypes)
        assert len(data) == num_dtypes, "Didn't get count dtypes %d" % num_dtypes
        if self.nbits == 2:
            data2 = np.zeros(num_elements*samps_per_element, dtype=np.int8)

            print('samp', num_samples, 'nelem', num_elements, 'ndtypes', num_dtypes, len(data), len(data2))

            data2[0::4] = (data & 0b00000011)*2 - 3
            data2[1::4] = (data & 0b00001100)/2 - 3
            data2[2::4] = (data & 0b00110000)/4 - 3
            data2[3::4] = (data & 0b11001100)/8 - 3

            data = data2
            
        data.shape = (num_samples, self.nchans)

        return data
    

    def __getitem__(self, slice_list):
        return self.get_data(slice_list, 0, 0)
        
    def print_header(self):
        for k,v, in self.header.items():
            if isinstance(v, float):
                print('{}:{:0.15f}'.format(k,v))
            else:
                print(k, ':', v)
        

def _main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_usage('%prog [options] args')
    parser.set_description('Script description')
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose [Default %default]')
    parser.set_defaults(verbose=False)
    (values, args) = parser.parse_args()
    global _verbose
    _verbose = bool(values.verbose)
    
    for filename in args:
        fin = SigprocFile(filename)
        fin.print_header()
        print('RADEC DEG', fin.src_raj_deg, fin.src_dej_deg)
        print("Header size", fin.header_size_bytes)
        print("Data size", fin.data_size_bytes)
        print("Number of elements", fin.file_size_elements)
        print('Duration (seconds)', fin.observation_duration)
    

if __name__ == '__main__':
    _main()
