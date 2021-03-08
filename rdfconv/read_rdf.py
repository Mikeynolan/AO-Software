#!/usr/bin/env python2
# Reads RDF-format radar data files
# Using Python 2
# Copied and modified from NAIC
#   /home/smarshal/Python/Read_RDF/read_rdf_v7b.py
# Relevant IDL scripts:
#   /home/cmagri/idl/io.pro (write*, rdfHeaderOK, readrdf)
#   /home/nolan/idl/cwcross.pro
#   And figure out what scripts get called to make RDFs when we are taking data
#     For instance, see log from 2018-03-22 survey night
#   IDL documentation:
#     http://www.harrisgeospatial.com/docs/using_idl_home.html
#     http://www.idlcoyote.com/
# By Sean Marshall
# Last modified on Thursday, May 17, 2018

import struct # file:///usr/share/doc/python-docs-2.7.5/html/library/struct.html
#import sys
import numpy
import matplotlib.pyplot as plt
# file:///usr/share/doc/python-docs-2.7.5/html/faq/programming.html#what-are-the-best-practices-for-using-import-in-a-module

interactive_mode = plt.isinteractive() # If running this script over a remote terminal, this should be False!
#print "interactive_mode :", interactive_mode
# Otherwise, plots won't work (over a remote terminal)
if not interactive_mode:
    import matplotlib
    matplotlib.use('Agg') # Allows plotting over remote terminal
    # To clarify: This will allow one to save all plots to image files,
    #             rather than displaying them in a new window or cell.
    # https://stackoverflow.com/questions/41814254/qxcbconnection-error-when-python-matplotlib-run-as-cron-job
    # file:///usr/share/doc/python-matplotlib-doc-1.2.0/api/matplotlib_configuration_api.html#matplotlib.use

# How an RDF file is organized:
# Header: format, name, value, comment (all in plain text, separated by whitespace)
#   For instance: "i       Height    500 # Number of range rows in image"
#   Format is always one lowercase character?
#   Casing of name may vary, e.g. could have "height" or "Height"
#   Comment is optional; if present, it begins with '#'
#   Lines/fields are separated by one byte, 0x0A
#   Required fields (I think):
#     For DD: type, height, width, size, machine, format
#     For CW: type, ndata, height, width, size, machine, nchan, ntags, format
# Header ends with three bytes: 0A 2E 0A (in hexadecimal; 10 46 10 in decimal)
# Data: binary values
#   Always stored as big-endian? Or just usually big-endian?
#   Typically (or always?) four-byte floating point
#   Each CW spectrum is typically followed by tags
#     These tags are stored as binary data, immediately after each channel's spectrum
#     But delay-Doppler images cannot be followed by tags?
#   For dual-polarization CW spectra, order is typically:
#     [OC spectrum 1, OC tags 1, SC spectrum 1, SC tags 1; OC spectrum 2, ...]
#   Chris Magri's io.pro (in comments near the beginning of function readrdf)
#     says that other orders are allowed. For instance:
#     [OC spectrum 1, OC tags 1, OC spectrum 2, OC tags 2; SC spectrum 1, ...]
#     or [OC 1, SC 1; SC 2, OC 2]
#   There are no markers (bytes) separating back-to-back spectra and tags
# There may (optionally) be one byte, 0x0A, separating the data (or tags) and the footer
#   Looks like images have that separator byte but CW spectra do not
# Footer (tail): format, name, value, comment (all in plain text, separated by whitespace)
#   Just like the header
#   Lines/fields are separated by one byte, 0x0A
#   For CW spectra, footer includes lines of text describing each binary tag
# Frame ends with three bytes: 0A 2E 0A
# Then the next frame's header begins (if there is another frame)

null_value = -1.0e20
#c_light = 299792458. # Speed of light, in meters per second
font_size = 16
spacer = "    " # Default indentation for some printed output

# Define constants and units (needed for describe_rdf_dict):
kB = 1.38064852e-23 # Boltzmann's constant: 1.380 648 52(79) * 10^-23 J/K
# From https://physics.nist.gov/cgi-bin/cuu/Value?k|search_for=all!
km = 1.0e3 # 1 kilometer is 1000 m
#kW = 1.0e3 # 1 kilowatt is 1000 W
Jy = 1.0e-26 # 1 Jansky is 10^-26 W/(m^2*Hz)
KJy = 2.0*kB/Jy
# 1 K/Jy is equivalent to an effective area of 2761. m^2 = 2*k_B*(1 K)/(1 Jy)
#   http://www.cv.nrao.edu/course/astr534/Equations.html
gain_unit_rdf = 1.0e12 # Conversion factor for dimensionless gain in RDF tags
# In an email from Chris (June 6, 2017):
#  "The "gain" tag value is (dimensionless TX gain)*(dimensionless RX gain)*(1e-12)
#   Put another way, it's 1e-12 times the product that shows up in the radar equation.
#   (I inherited this obscure usage so I take no responsibility for it.)
#   Or, if you prefer, for Arecibo data it's
#   (2.1846 * TX sensitivity in K/Jy)*(2.1846 * RX sensitivity in K/Jy)
#   where the "2.1846" includes various powers of ten and factors of pi and
#   Boltzmann's constant and the observing wavelength (so you'd need a different
#   number for Goldstone data). Your value of 440 corresponds to 9.6 K/Jy."

def ainfo(inp, spacing=""):
    """Prints information on input NumPy array"""
    # This will be used by the function describe_rdf_dict (defined below)
    print spacing + "Array type:", type(inp)
    if inp.ndim == 0:
        temp = type(numpy.atleast_1d(inp)[0])
        # Handling zero-rank arrays:
        #   https://stackoverflow.com/questions/4565749/how-detect-length-of-a-numpy-array-with-only-one-element
        #   https://web.archive.org/web/20130401215434/http://projects.scipy.org/numpy/wiki/ZeroRankArray
        #   https://stackoverflow.com/questions/773030/why-are-0d-arrays-in-numpy-not-considered-scalar
        #   https://stackoverflow.com/questions/16805987/python-numpy-ndarray-subclasses-and-zero-rank-arrays
        #   https://stackoverflow.com/questions/17414098/numpy-zero-rank-array-indexing-broadcasting
    else:
        temp = type(inp.flat[0])
        # Could use inp.dtype instead of type(inp.flat[0])
        # But that would show something like "float64" instead of "<type 'numpy.float64'>"
    print spacing + "Element type:", temp, \
          "({0:,d} bytes each)".format(inp.itemsize)
    print spacing + "Size:", inp.size, "({0:,d} bytes total)".format(inp.nbytes)
    print spacing + "Shape:", inp.shape
    #return inp.shape
    pass

# Function is_text(chunk) used to be here

def read_file(file_path):
    """Reads a file of the radar data format (RDF); returns a bytearray"""
    #"""Reads a file of the radar data format (RDF); returns NumPy uint8 array"""
    inp = numpy.fromfile(file_path, dtype=numpy.uint8)
    #return inp
    return bytearray(inp)

new_line = bytearray([10]) # Indicates new line in RDF header or footer
new_section = bytearray([10, 46, 10]) # Indicates end of RDF header or footer
new_section_len = len(new_section)

rdf_search_strings = ["s                 type   float", \
                      "s          Type   Float", \
                      "s       Type      Float", \
                      "c       Type      Float"]
# Strings that can indicate the start of a frame (section) of an RDF
#   From various RDF files that I encountered during testing

def find_rdf_frames(inp):
    """Determines whether the input RDF data block (bytearray) contains more
    than one frame (spectrum or image);
    Returns addresses of each section's (frame's) beginning and end"""
    # inp should be a bytearray
    unknown_search_string = True
    # There are several possibilities to check (in rdf_search_strings), but this
    #   function assumes that any input RDF uses just one of the possible
    #   strings (rather than using different strings for different frames).
    ics = 0 # Counter for entries of rdf_search_strings
    #while unknown_search_string and ():
    for ics in range(len(rdf_search_strings)):
        # Try each entry in rdf_search_strings; see if it is in the input array
        # TODO: Do this search using regular expressions instead of fixed strings?
        #       How does Chris Magri handle this in his io.pro ?
        N_frames = inp.count(rdf_search_strings[ics])
        if N_frames > 0:
            #search_string = rdf_search_strings[ics]
            unknown_search_string = False
            if N_frames == 1:
                # If there's only one frame, just return beginning and ending addresses
                return [ 0, len(inp) ]
            else:
                # If there's more than one frame
                indices = [ -1 ] # Will append to this list to track where the search string appears
                continue_search = True
                while continue_search:
                    try:
                        indices.append(inp.index(rdf_search_strings[ics], \
                                                 indices[-1]+1))
                    except ValueError:
                        continue_search = False
                indices.append(len(inp)) # Ending address for the last frame
                return indices[1:]
    #if N_frames <= 0:
    if unknown_search_string:
        raise ValueError("Cannot find any of rdf_search_strings in input data")
    # The function should never get past this point
    print "Something weird has happened, in function read_rdf.find_rdf_frames"
    pass

def parse_header_line(inp, verbose=False):
    """Parses one line from an RDF's header or footer
    Returns tuple of strings, (format, name, value, comment)
    Even if the value is numeric, it is returned as a string."""
    if verbose:
        print spacer + "parse_header_line:", inp
    # TODO: Handle vectors of values, e.g. "v       stat[5] 14 12 13 14 0"
    comment_count = inp.count('#')
    if comment_count > 1:
        raise ValueError("Input string has multiple instances of # character")
    comment_ind = inp.find('#') # -1 if there's no instance of '#'
    inp_split = inp.split() # Split at whitespace
    format_char = inp_split[0]
    if len(format_char) != 1:
        raise ValueError("Format substring should be a single character")
    #if comment_ind == 0:
    if format_char == '#':
        # If inp has ONLY a comment, e.g. "#       writerdf  Revision: 1.10 "
        name_str = "" # No name or value is present
        value_str = "" # So use the empty string for those
        comment_str = inp[comment_ind:]
        # Keep comment's leading '#' character
        # Also keep comment's trailing whitespace, if there is any (no strip)
    else:
        # Normal header/footer line, with format, name, value, and optional comment
        #   e.g. "i       Width     200 # Number of Doppler columns in image"
        name_str = inp_split[1] # A field's name should only be one word
        # The value can include spaces, so it must be handled more carefully
        if comment_count == 0:
            # value is everything after name
            value_str = inp[inp.index(inp_split[2]):].strip()
            comment_str = "" # If no comment is present; use the empty string
        else:
            # value is everything after name and before the comment
            value_str = inp[inp.index(inp_split[2]):comment_ind].strip()
            comment_str = inp[comment_ind:]
            # Keep comment's leading '#' character
            # Also keep comment's trailing whitespace, if there is any (no strip)
    #return (format_char, name_str, value_str, comment_str)
    return format_char, name_str, value_str, comment_str

def parse_rdf(inp, file_path="Source file not specified", frame=-32768):
    """Parses the contents of a file of the radar data format (RDF)
    Input should be a bytearray
    Returns a list of dictionaries (one dict for each frame in the RDF)"""
    inp_len = len(inp)
    outp = [('source_file', file_path), ('source_file_frame', frame)]
    comments = [ ] # Will store comment-only lines from header and footer
    height = 0 # Number of delay rows, or number of CW spectra
    width = 0 # Number of frequency bins per row
    ndata = -32768 # Number of non-tag data points in each image or spectrum
    bytes_per_point = 0 # Size of each data value
    ntags = 0 # Number of tags that are stored as binary data (NOT text)
    tag_names = [ ] # Will store names of tags
    tag_indices = [ ] # Will store indices of tags
    current_frame_number = 0 # Zero-based index
    address_start = 0 # Start of file
    address_end   = inp.index(new_line)
    section_end = inp.index(new_section)
    addresses = [("Frame {0:d} header".format(current_frame_number), \
                  address_start)]
    data_start = section_end + new_section_len
    # index raises a ValueError if it can't find the search substring
    while address_start < inp_len - new_section_len:
        format_char, name_str, value_str, comment_str = \
                         parse_header_line(str(inp[address_start:address_end]))
        if format_char == 'c' or format_char == 's':
            if name_str.lower() == "type":
                if value_str.lower() != "float":
                    raise ValueError("Data type is " + value_str + \
                                     "; should be float")
            outp.append((name_str, [value_str, format_char, comment_str]))
        elif format_char == '#':
            comments.append(comment_str)
        elif format_char == 'i':
            outp.append((name_str, [int(value_str), format_char, comment_str]))
            if name_str.lower() == "ndata":
                ndata = int(value_str)
            elif name_str.lower() == "height":
                height = int(value_str)
            elif name_str.lower() == "width":
                width = int(value_str)
            elif name_str.lower() == "size":
                bytes_per_point = int(value_str)
            elif name_str.lower() == "ntags":
                ntags = int(value_str)
        elif format_char == 'f' or format_char == 'd':
            outp.append((name_str, \
                         [float(value_str), format_char, comment_str]))
        # TODO: CONTINUE HERE. Is format_char always lowercase?
        elif format_char == 't':
        elif format_char == 'v':
        else:
            print "Error parsing line:"
            print str(inp[address_start:address_end])
            raise ValueError("Unrecognized format")
        address_start = address_end + 1
        address_end  = inp.index(new_line, address_start, section_end)
        pass
    outp.append(('comments', comments))
    raise NotImplementedError("Function under construction")

# TODO: Should parsing each frame of an RDF be in one function, and then
#         finding frames would be in another function?
# Or have parse_rdf return something indicating whether it reached the end of the file

def parse_rdf_OLD(inp, file_path="Source file not specified", frame=-32768):
    """Parses the contents of a file of the radar data format (RDF), returning dict"""
    # inp should be a NumPy array of uint8 values
    # Find where inp has bytes with value 0x0a (decimal 10)
    # Since those may indicate separations between lines of text
    ind010 = numpy.where(inp == 10)[0] # Indices separating possible chunks of text
    # Now process those chunks
    chunk = str(bytearray(inp[0:ind010[0]]))
    chunk_spl = chunk.split()
    temp = chunk.find(chunk_spl[2])
    outp = [('source_file', file_path), ('source_file_frame', frame), \
            (chunk_spl[1], chunk[temp:])]
    text = [ chunk ] # Also keep a separate list with all text
    # First chunk should be text (NOT binary data)
    # Since some numeric tags include comments that may be needed later
    height = 0 # Number of delay rows, or number of CW spectra
    width = 0 # Number of frequency bins per row
    bytes_per_point = 0 # Size of each data value
    ntags = 0 # Number of tags that are stored as binary data (NOT text)
    tag_names = [ ] # Will store names of tags
    tag_indices = [ ] # Will store indices of tags
    ndata = -32768 # Number of non-tag data points in each image or spectrum
    #data_start = [ ] # Will store the starting positions of each block of data
    ici = 0 # Counter for all chunks
    icd = 0 # Counter for data (non-text) chunks that have been parsed
    for ici in range(ind010.size - 1):
        # The final chunk should be just a tab, so skip it
        chunk = bytearray(inp[ind010[ici]:ind010[ici+1]])
        if is_text(chunk[1:]):
            chunk = str(chunk[1:]) # Skip the first character (newline)
            text.append(chunk)
            #print chunk # For debugging
            chunk_spl = chunk.split()
            if chunk[0] == 'c' or chunk[0] == 's':
                temp = chunk.find(chunk_spl[2])
                outp.append((chunk_spl[1], chunk[temp:]))
            elif chunk[0] == 'v':
                # Vector of values, e.g. "v       stat[5] 14 12 13 14 0"
                #print chunk # For debugging
                ind_temp = [chunk.find(chunk_spl[1]), \
                            chunk.find('['), chunk.find(']')]
                #print chunk[ind_temp[0]:ind_temp[1]], chunk[ind_temp[2]+1:]
                #outp.append((chunk[ind_temp[0]:ind_temp[1]], chunk[ind_temp[2]+1:]))
                temp = chunk[ind_temp[2]+1:].split()
                v_values = [ ]
                icv = 0 # Counter for vector values
                for icv in range(len(temp)):
                    v_values.append(float(temp[icv]))
                #print v_values
                if int(chunk[ind_temp[1]+1:ind_temp[2]]) != len(v_values):
                    print "WARNING: Unexpected length for vector " + \
                          chunk[ind_temp[0]:ind_temp[1]]
                outp.append((chunk[ind_temp[0]:ind_temp[1]], \
                             numpy.array(v_values)))
            elif chunk[0] == 't':
                # Name of numeric tag
                # Tag values are stored as binary data after each "row" of data
                tag_names.append(chunk_spl[1])
                tag_indices.append(int(chunk_spl[2]))
                #print "  New tag: " + chunk_spl[1] # For debugging
            elif chunk[0] == 'd':
                outp.append((chunk_spl[1], float(chunk_spl[2])))
            elif chunk[0] == 'i':
                outp.append((chunk_spl[1], int(chunk_spl[2])))
                if chunk_spl[1].lower() == "height":
                    height = int(chunk_spl[2])
                elif chunk_spl[1].lower() == "width":
                    width = int(chunk_spl[2])
                elif chunk_spl[1].lower() == "size":
                    bytes_per_point = int(chunk_spl[2])
                elif chunk_spl[1].lower() == "ntags":
                    ntags = int(chunk_spl[2])
                elif chunk_spl[1].lower() == "ndata":
                    ndata = int(chunk_spl[2])
            elif chunk[0] == 'f':
                outp.append((chunk_spl[1], float(chunk_spl[2])))
            elif chunk[0] == '#':
                pass
            elif chunk[0] == ' ':
                pass
            else:
                print "Error in chunk:"
                print chunk
                raise ValueError("Unrecognized tag format")
        else:
            icd += 1
            #print icd, ind010[ici], ind010[ici+1]
            #print inp[ind010[ici]:ind010[ici]+4]
            if icd == 2:
                row_start = ind010[ici] + 1
                if ntags > 0:
                    # Need to handle initial tag separately
                    # Since it is not preceded by 0x0a
                    if len(tag_names) > 0:
                        raise ValueError("Name of tag 0 was set before data" + \
                                         " block")
                    if height > 0 and width > 0:
                        ind_tag0 = row_start + height*width*bytes_per_point
                        ind_tag1 = ind010[numpy.min( \
                                          numpy.where(ind010 > ind_tag0)[0])]
                        #print "Addresses for start and end of tag 0:", ind_tag0, ind_tag1
                        chunk = bytearray(inp[ind_tag0:ind_tag1])
                        if is_text(chunk):
                            chunk = str(chunk) # No leading newline to skip
                            text.append(chunk)
                            #print chunk # For debugging
                            chunk_spl = chunk.split()
                            if chunk[0] == 't':
                                tag_names.append(chunk_spl[1])
                                tag_indices.append(int(chunk_spl[2]))
                                #print "  New tag: " + chunk_spl[1] # For debugging
                        else:
                            raise ValueError("No text at expected address of" + \
                                             " first tag")
                    else:
                        raise ValueError("Data block began before height" + \
                                         " and width were set")
            #else:
            #    if len(chunk) > 2:
            #        data_end = ind010[ici] + 1
    if height <= 0:
        raise ValueError("Problem with height")
    if width <= 0:
        raise ValueError("Problem with width")
    if bytes_per_point <= 0:
        raise ValueError("Problem with bytes_per_point")
    if ndata < 0:
        #print "WARNING: ndata had been set to", ndata
        ndata = width - ntags
    if ndata + ntags != width:
        raise ValueError("Problem with ndata")
    #data_start.append(data_start[0] + width*bytes_per_point)
    #data_start = numpy.array(data_start)
    row_start = row_start + numpy.arange(height)*width*bytes_per_point
    row_end = row_start + ndata*bytes_per_point
    #data_end = data_start + height*ndata*bytes_per_point
    #print data_start, data_end, data_end - data_start
    #print inp[data_start:data_start+4]
    #print inp[data_end-4:data_end]
    #outp.append(('data', inp[data_start:data_end]))
    #N_pts = height*width
    data = null_value*numpy.ones((height, ndata))
    tag_values = null_value*numpy.ones((height, ntags))
    #      Initialize these arrays by filling them with an invalid value
    ici = 0 # Counter for rows of D-D images, or for which spectrum in a series
    #ick = data_start - bytes_per_point # Counter for address of inp
    for ici in range(height):
        icj = 0 # Counter for frequency (in D-D pixels or frequency bins), then for tags
        ick = row_start[ici] - bytes_per_point
        for icj in range(width):
            ick += bytes_per_point
            # TODO: Do this without needing one step of the loop for each point
            if ick < row_end[ici]:
                # If it's still in the data block
                data[ici,icj] = struct.unpack('>f', \
                                              inp[ick:ick+bytes_per_point])[0]
                # Previous line assumes float values are stored in big-endian format!
            else:
                # If it's in the tags
                tag_values[ici,icj-ndata] = struct.unpack('>f', \
                                                inp[ick:ick+bytes_per_point])[0]
    #print tag_names # For debugging
    #print tag_indices
    #ainfo(tag_values)
    #print tag_values
    icj = 0 # Counter for tags
    if ntags > 0:
        if height <= 1:
            #tag_values = tag_values.flatten()
            for icj in range(ntags):
                outp.append((tag_names[icj], tag_values[0,tag_indices[icj]]))
        else:
            for icj in range(ntags):
                #print tag_names[icj], tag_indices[icj], tag_values[:,tag_indices[icj]]
                outp.append((tag_names[icj], tag_values[:,tag_indices[icj]]))
    temp = dict(outp)
    if ("dfreq" in temp.viewkeys()) and ("xjcen" in temp.viewkeys()):
        if height <= 1:
            frequencies = temp['dfreq']*(numpy.arange(ndata) - temp['xjcen'])
            # xjcen should be zero-based, since IDL (like Python) uses
            #   zero-based indices. In my earlier version of this script
            #   (renamed to read_rdf_v6_OLD1_BAD.py), I decremented
            #   temp['xjcen']. That was a mistake; now it's correct.
            # This index correction is new to v6! I missed it in v5 and earlier!
        #elif height == 2:
        elif height >= 2:
            frequencies = temp['dfreq'][0]*(numpy.arange(ndata) - \
                                            temp['xjcen'][0])
            if temp['xjcen'][0] != temp['xjcen'][1]:
                #print "WARNING: xjcen differs for the first two channels"
                print "WARNING: xjcen differs for the first two spectra"
            if temp['dfreq'][0] != temp['dfreq'][1]:
                print "WARNING: dfreq differs for the first two spectra"
        #else:
        #    # I don't think I have ever seen an RDF like this
        #    raise ValueError("height = {0:d}".format(height))
        outp.append(('frequencies', frequencies))
    threshold = 10000.0
    temp = numpy.min(data)
    if temp < -threshold:
        print "WARNING: Minimum data value is {0:+.3e}".format(temp)
        print spacer + \
              "Data may not be properly normalized, but attempting to continue"
    temp = numpy.max(data)
    if temp > threshold:
        print "WARNING: Maximum data value is {0:+.3e}".format(temp)
        print spacer + \
              "Data may not be properly normalized, but attempting to continue"
    outp.append(('data', data))
    outp.append(('text', text))
    return dict(outp)

def load_rdf(file_path, verbose=True):
    """Function to read an RDF file, find its frame boundaries, and
    return a dict (for one frame) or list of dicts (for multiple frames)"""
    # This is intended to be a simple way of reading an RDF file and
    #   getting it into an easily usable format, with just one function call.
    # How to call each function:
    # read_file(file_path_string)
    # find_rdf_frames(bytearray)
    # parse_rdf(inp_NumPy_uint8_array, file_path="Source file not specified", frame=-32768)
    if verbose:
        print "Reading data from file", file_path
    inp = read_file(file_path)
    # inp is a numpy.ndarray of numpy.uint8 values
    if verbose:
        print "Info for binary data from loaded file:"
        ainfo(inp, spacer)
    sect_addr = find_rdf_frames(bytearray(inp))
    if verbose:
        print "Addresses of frame boundaries:"
        #print numpy.array(sect_addr), type(sect_addr), len(sect_addr)
        print sect_addr, type(sect_addr), len(sect_addr)
    outp = [ ] # Will append each frame's dict to this list
    # Each element of this list will be a parse_rdf dict from a different frame
    # But, both polarizations (if present) from a given frame go into the same dict
    icf = 0 # Counter for frames (within the loaded RDF file)
    for icf in range(len(sect_addr) - 1):
        outp.append(parse_rdf(inp[sect_addr[icf]:sect_addr[icf+1]], \
                              file_path, icf))
    if len(sect_addr) == 2:
        # If sect_addr has two items, the loaded RDF only has data from one frame
        # Then sect_addr[0] should be 0 (start of file)
        #  and sect_addr[1] should be the file size in bytes (meaning end of file)
        #outp = parse_rdf(inp[sect_addr[0]:sect_addr[1]], file_path, 0)
        if verbose:
            print "The loaded RDF has data from one frame."
            print "This output returned by this function will be a list."
            print "The list will have one item, a dictionary (dict)."
            print "Most of the dictionary's keys correspond to the RDF tags."
            if "frequencies" in outp[0].viewkeys():
                print "Also, output[0]['frequencies'] has the frequencies,"
                print "  and output[0]['data'] has the data (CW spectra)."
            else:
                print "Also, output[0]['data'] has the data (D-D image)."
    elif len(sect_addr) > 2:
        # If sect_addr has more than two items, the loaded RDF has data from multiple frames
        # Then sect_addr[0] should be 0 (start of file, and of first frame)
        #      sect_addr[1] should be the start of the second frame
        #      sect_addr[2] should be the start of the third frame
        #      ... and so on ...
        #      sect_addr[-1] should be the file size in bytes (meaning end of file)
        #outp = [ ] # Will append each frame's dict to this list
        if verbose:
            #print type(outp), len(outp) # outp is a list
            #print type(outp[0]), len(outp[0].keys()) # Each element of outp is a dictionary
            # outp[0] has the first frame from the RDF file
            # outp[1] has the second frame from the RDF file
            # outp[2] has the third frame from the RDF file, and so on...
            print "The loaded RDF has data from", len(sect_addr)-1, "frames."
            print "This output returned by this function will be a list."
            print "Each item in that list will be a dictionary (dict) from" + \
                  " a different frame."
            print "    output[0] is the dictionary for the first frame,"
            print "    output[1] is the dictionary for the second frame, etc."
            print "Most of a dictionary's keys correspond to the RDF tags."
            if "frequencies" in outp[0].viewkeys():
                print "Also, output[#]['frequencies'] has the frequencies,"
                print "  and output[#]['data'] has the data (CW spectra)."
            else:
                print "Also, output[#]['data'] has the data (D-D image)."
            #print sorted(list(outp[0].viewkeys()))
    else:
        # I don't think this should ever happen, since find_rdf_frames should
        #   raise an error first, if it can't find any frames/sections.
        raise ValueError("Problem finding frame boundaries")
    return outp

# For each entry (sub-list) in cw_lists: First element is RDF field name,
#   second element is description (before value),
#   third element is text for after value(s) (usually specifying units),
# TODO: Add final element that gives print format specification?
# I did something similar in my ~/Python/horizcheck.py
cw_lists = [ ['source_file', "Source file:", ""] ]
cw_lists.append(['source_file_frame', \
                 "  Frame (section) of source file:", "(zero-based index)"])
cw_lists.append(['format', "  Format:", ""])
cw_lists.append(['type', "  Data type:", ""])
cw_lists.append(['size', "  Data point size:", "bytes (per point)"])
cw_lists.append(['machine', "  Machine:", ""])
cw_lists.append(['target', "Target name:", ""])
cw_lists.append(['xmit_sta', "Transmitter station:", ""])
cw_lists.append(['calmean', "Mid-receive time:", ""])
cw_lists.append(['iyy',    "  Year:  ", ""])
cw_lists.append(['imm',    "  Month: ", ""])
cw_lists.append(['idd',    "  Day:   ", ""])
cw_lists.append(['rchour', "  Hour:  ", ""])
cw_lists.append(['rcmin',  "  Minute:", ""])
cw_lists.append(['rcsec',  "  Second:", ""])
cw_lists.append(['rcnsec', "  Nanos.:", ""])
cw_lists.append(['timezone', "  Time zone:", ""])
cw_lists.append(['tzcorr', "  Time zone correction:", ""])
cw_lists.append(['rcsta', \
                 "  Receive start time:", "seconds after UT midnight"])
cw_lists.append(['rcend', \
                 "  Receive end time:  ", "seconds after UT midnight"])
cw_lists.append(['jd0',     "  JD0:           ", ""])
cw_lists.append(['jdstart', "  Start JD:      ", ""])
cw_lists.append(['jdmean',  "  Mid-receive JD:", ""])
cw_lists.append(['jdend',   "  Stop JD:       ", ""])
cw_lists.append(['elev', "Target's elevation:", "degrees"])
cw_lists.append(['azim', "Target's azimuth:  ", "degrees"])
#cw_lists.append(['zepch', "Epoch:", ""]) # This may be something else
cw_lists.append(['ramin',  "Minimum R.A.:", "hours"])
cw_lists.append(['ramean', "Mean R.A.:   ", "hours"])
cw_lists.append(['ramax',  "Maximum R.A.:", "hours"])
cw_lists.append(['decmin',  "Minimum declination:", "degrees"])
cw_lists.append(['decmean', "Mean declination:   ", "degrees"])
cw_lists.append(['decmax',  "Maximum declination:", "degrees"])
cw_lists.append(['distmin',  "Minimum distance:", "au"])
cw_lists.append(['distmean', "Mean distance:   ", "au"])
cw_lists.append(['distmax',  "Maximum distance:", "au"])
cw_lists.append(['rttim', "Round-trip light time:", "seconds"])
cw_lists.append(['tau',   "Receive duration:     ", "seconds"])
cw_lists.append(['lambda', "Radar wavelength:", "meters"])
cw_lists.append(['gain', "Gain:", "(RDF units; two-way)"])
cw_lists.append(['trpwr', "Transmitter power:", "kW"])
cw_lists.append(['tsys', "System temperature:", "K"])
cw_lists.append(['sdev', "Noise-equivalent cross section:", "km^2"])
# sdev usually won't print clearly (with my default NumPy print options)
#   Because its values are typically too small
cw_lists.append(['nchan',  "Number of channels:", ""])
cw_lists.append(['nspec',  "Number of spectra: ", ""])
cw_lists.append(['height', "Height:            ", ""])
#                 Number of delay rows, or number of CW spectra?
# The values for 'height', 'nchan', and 'nspec' should all be equal
cw_lists.append(['width', "Width:                        ", ""])
# Number of frequency bins per row (plus number of tags)?
cw_lists.append(['ndata', "Number of non-tag data points:", ""])
#cw_lists.append(['ntags', "Number of tags that are stored as binary data (NOT text):", ""])
cw_lists.append(['ntags', "Number of (binary data) tags: ", ""])
# The value of width should be equal to (ndata + ntags)
cw_lists.append(['lfft', "FFT length:", ""])
cw_lists.append(['xjcen', "Ephemeris COM channel:", "(zero-based index)"])
# xjcen should be equal for the two channels!
cw_lists.append(['doppl', "Doppler offset:", "Hz"])
cw_lists.append(['dfreq', "Frequency resolution:", "Hz"])
#cw_lists.append(['kpts', "Data points kept", ""]) # Unsure about this one
cw_lists.append(['nffts', "Number of looks:", ""])
# nffts is the number of looks, according to Mike's cwcross.pro
# But, based on Chris and Mike's IDL scripts, there can also (possibly?) be another tag called nlooks
#cw_lists.append([, "", ""])
cw_lists.append(['frequencies', "Frequencies:", "Hz"])
cw_lists.append(['data', "Data values:", "(SNR; unitless)"])
cw_lists.append(['text', "Text from parsing RDF:", ""])
N_cw_fields_arrays = 3 # Number of entries that are (large) arrays or lists
# Those will NOT be printed in their entirety
# Fields from the test RDF (1999 JD6 from Arecibo on 2015-07-29) that I don't recognize:
#[, , 'color', 'crerr', 'cross', , , ,\
# , , 'diameter', , , , , \
# , , 'freq1', , 'frstep', , , \
# , 'igw', , 'irun', 'itar', , 'jcp', , , \
# , , 'jgroup', 'jsnr1', 'jsnr2', 'kpts', , ,\
# , , , , 'nfreq', , , 'obs', \
# 'period', 'phase', 'phase0', 'posfr', , , , ,\
# , , , , , 'rmsc', 'rmsm', , \
# , , , , , , \
# , , , , , , 'util', , \
# , , ]
# cross should be the (automatically computed) cross section?
#   And crerr should be its uncertainty?
# jcp indicates polarization. Normally it's 0 for OC and 1 for SC.
#   But, according to Chris's io.pro, tkplay uses 1 for OC and 2 for SC
# jsnr? are indices bounding signal, for cross section calculations?
# kpts is number of data points acquired? (12,500 samples per second)*(Rx duration)
# posfr is related to frequency limits in Chris's plots
# diameter and period should be filled in automatically?
#   phase0 and phase are derived from the period and the observation time? And jd0?
#cw_names = list(item[0] for item in cw_lists) # Names of fields in the dict from parse_rdf
#cw_descriptions = list(item[1] for item in cw_lists) # Description strings
cw_names = [ ]
cw_descriptions = [ ]
cw_units = [ ]
for item in cw_lists:
    # Fill in the lists that will be used for printed output
    cw_names.append(item[0])
    cw_descriptions.append(item[1])
    cw_units.append(item[2])
    # TODO: Use built-in function zip instead of this loop? (but that returns tuples?)
    #       file:///usr/share/doc/python-docs-2.7.5/html/library/functions.html#zip
    #       file:///usr/share/doc/python-docs-2.7.5/html/faq/programming.html#my-program-is-too-slow-how-do-i-speed-it-up

def describe_rdf_cw_dict(inp_dict, spacing="  "):
    """Given a dict from parse_rdf, prints many fields' values and descriptions
    This was written and tested with CW spectra.
    It won't work well with delay-Doppler images (due to different field names)."""
    inp_keys = sorted(list(inp_dict.viewkeys()))
    inp_keys_lowercase = [ ] # Will store lowercase versions of all inp_keys items
    # To enable case-insensitive searches
    for item in inp_keys:
        inp_keys_lowercase.append(item.lower())
        # Elements of inp_keys and inp_keys_lowercase are in the same order
    print "Known tags:"
    N_fields = len(cw_names)# - 1 # Since the last item is a placeholder
    ici = 0 # Counter
    for ici in range(N_fields - N_cw_fields_arrays):
        #if cw_names[ici] in inp_keys:
        if cw_names[ici] in inp_keys_lowercase:
            # If this key is in the input dict
            icj = inp_keys_lowercase.index(cw_names[ici])
            # Searching cw_names for each key seems inefficient
            # But it's easy, and there aren't very many elements in cw_names
            print spacing + cw_descriptions[ici], inp_dict[inp_keys[icj]], \
                  cw_units[ici] + " ['" + inp_keys[icj] + "']"
            #      cw_units[ici] + " (" + cw_names[ici] + ")"
            if cw_names[ici] == 'gain':
                # For gain, also print related quantities
                #temp = numpy.sqrt(inp_dict['gain']*gain_unit_rdf) # Dimensionless gain
                temp = numpy.sqrt(inp_dict[inp_keys[icj]]*gain_unit_rdf)
                print spacing + spacer, temp, "(dimensionless; one-way)"
                print spacing + spacer, 10.0*numpy.log10(temp), "dB (one-way)"
                # G = 4*pi*A_eff/lambda^2, so A_eff = G*lambda^2/(4*pi)
                #temp *= inp_dict['lambda']**2/(4.0*numpy.pi) # Effective area, in m^2
                icj = inp_keys_lowercase.index('lambda')
                # list.index() should raise a ValueError if there is no matching item
                temp *= inp_dict[inp_keys[icj]]**2/(4.0*numpy.pi)
                print spacing + spacer, temp/KJy, "K/Jy"
                print spacing + "  Effective area:", temp, "m^2"
                # A = pi*R^2 = (pi/4)*D^2; D^2 = 4*A/pi
                print spacing + "  Effective diameter:", \
                      numpy.sqrt(4.0*temp/numpy.pi), "m"
            if cw_names[ici] == 'sdev':
                # For sdev, also print it in m^2
                #print spacing + 30*" " + "=", inp_dict['sdev']*km*km, "m^2"
                print spacing + 30*" " + "=", \
                      inp_dict[inp_keys[icj]]*km*km, "m^2"
        else:
            print "WARNING: Expected key " + cw_names[ici] + \
                  " was not in the input dict"
    for ici in range(N_fields - N_cw_fields_arrays, N_fields):
        if cw_names[ici] in inp_keys_lowercase:
            icj = inp_keys_lowercase.index(cw_names[ici])
            if cw_names[ici] == 'text':
                # Chunks of text that were found when parsing the RDF
                print spacing + cw_descriptions[ici], \
                      type(inp_dict[inp_keys[icj]]), \
                      len(inp_dict[inp_keys[icj]]), \
                      cw_units[ici] + " ['" + cw_names[ici] + "']"
            else:
                print spacing + cw_descriptions[ici], \
                      cw_units[ici] + " ['" + inp_keys[icj] + "']"
                ainfo(inp_dict[inp_keys[icj]], spacer)
        else:
            print "WARNING: Expected key " + cw_names[ici] + \
                  " was not in the input dict"
    print "Unknown tags:"
    for item in inp_keys:
        if item.lower() in cw_names:
            # If it is one of the known tags from above, skip it
            pass
        else:
            print spacing + item + ":", inp_dict[item]
    pass

def plot_rdf_img(inp, output_plot_fname=False):
    """Given a dictionary from read_rdf, displays data"""
    # http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.imshow
    # http://matplotlib.org/examples/color/colormaps_reference.html
    fig = plt.figure(figsize=(12, 10))
    plt.imshow(inp['data'], cmap='gray')
    plt.title(inp['source_file'].split('/')[-1], fontsize=font_size+8)
    plt.xlabel("Doppler column", fontsize=font_size+4)
    plt.ylabel("Delay row", fontsize=font_size+4)
    cbar = plt.colorbar()
    cbar.set_label("SNR", fontsize=font_size+4)
    cbar.ax.tick_params(labelsize=font_size)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size)
    fig.tight_layout()
    if output_plot_fname:
        # If given a second input argument, save the plot to that file name
        plt.savefig(output_plot_fname)
    else:
        plt.show()
    plt.close()
    #return output_plot_fname
    pass

def plot_rdf_cw(inp, maxf=20.0, output_plot_fname=False):
    """Given a dictionary from read_rdf, displays CW spectra"""
    fig = plt.figure(figsize=(12, 10))
    #N_freqs = inp['data'].shape[1] + 0
    #N_freqs = inp['ndata'] + 0
    #freqs = inp['dfreq'][0]*(numpy.arange(N_freqs) - inp['xjcen'][0])
    #if inp['xjcen'][0] != inp['xjcen'][1]:
    #    print "WARNING: xjcen differs for the two channels"
    multiple_spectra = False
    #if "nchan" in inp.viewkeys():
    #    if inp['nchan'] > 1:
    if "nspec" in inp.viewkeys():
        if inp['nspec'] > 1:
            #multiple_channels = True
            multiple_spectra = True
            plt.plot(inp['frequencies'], inp['data'][0,:], '.b')
            plt.plot(inp['frequencies'], inp['data'][1,:], '.r')
        else:
            plt.plot(inp['frequencies'].flatten(), inp['data'].flatten(), '.b')
        if inp['nspec'] > 2:
            print "WARNING: Only plotting the first two spectra" + \
                  " (out of {0:d})".format(inp['nspec'])
    else:
        plt.plot(inp['frequencies'].flatten(), inp['data'].flatten(), '.b')
    plt.title(inp['source_file'].split('/')[-1], fontsize=font_size+8)
    #plt.xlabel("Doppler column", fontsize=font_size)
    plt.xlabel("Doppler frequency (Hz)", fontsize=font_size)
    plt.ylabel("Signal-to-noise ratio", fontsize=font_size)
    if multiple_spectra:
        #plt.legend(("OC", "SC"), loc='best', framealpha=0.75, numpoints=1, \
        plt.legend(("0", "1"), loc='best', framealpha=0.75, numpoints=1, \
                   fontsize=font_size)
    if maxf > 0:
        plt.xlim((-maxf, maxf))
    plt.grid(True)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size)
    fig.tight_layout()
    if output_plot_fname:
        # If given a third input argument, save the plot to that file name
        plt.savefig(output_plot_fname)
    else:
        plt.show()
    plt.close()
    #return output_plot_fname
    pass

# Example usage: See read_rdf_v6a.py
#                and ../CW_cross_sections/cw_cross_v6a_RDF_test?.ipynb

# TODO: Go through Chris's io.pro, and follow his procedures
#         Continue looking at his readrdf (starting on line 2939)
#         Make sure that the output from my parse_rdf function includes a field
#           that lists the binary tag names (in their original order)
#           And do the same for footer?
#           And include a list (of tuples?) that specifies addresses in the original RDF file
#         NotImplementedError if the spectra are in a weird order
#           That is, something other than [OC, SC; SC, OC; ...]
#         Write a function analogous to his rdfHeaderOK
#         Follow his nomenclature
#           e.g. "tags" only refers to values stored as binary after the data (?)
#             Tags only occur for CW data, not for delay-Doppler images?
#             "tail" and "footer" mean plain-text format/labels/values at the end?
#       When loading delay-Doppler images, include fields (1-D arrays) for the
#         delay and Doppler values, analagous to 'frequencies' for CW
#       Write a function describe_rdf_dd_dict (for delay-Doppler images)
#         See writefitsi in /home/cmagri/idl/io.pro
#       Check whether this script still works for all previous test cases
#         Since I only tested this version with CW RDF files for (85989) 1999 JD6
#         Check the files that were loaded for ivar's /data/smarshall/Documents/FDL/Read_RDF/Read_RDF_v4.ipynb
#         And also NAIC /share/radar1/lbenner/phaethon/phaethon.dec08.runs.1Hz.rdf
#       Include the functions that calculate Phil Perillat's S-band narrow gain curves?
#         e.g. from /home/smarshal/Python/CW_cross_sections/cw_cross_v6a_RDF_test1.ipynb
#       Allow a verbose option, to print more debugging output?
#       Eventually, go through the various Magri/Nolan/??? scripts and figure out what the unknown RDF tags mean
#       Define a class for loaded RDF files?
#         https://en.wikipedia.org/wiki/Object-oriented_programming#External_links
#         file:///usr/share/doc/python-doc/html/tutorial/classes.html
#         file:///usr/share/doc/python-doc/html/reference/datamodel.html#new-style-and-classic-classes
#         file:///usr/share/doc/python-doc/html/reference/compound_stmts.html#class
