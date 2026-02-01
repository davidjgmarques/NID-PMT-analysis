/**
 * @file
 * @brief Lecroy waveform binary reader, standalone header file
 * @note Tested with binary template format 2_3
 */

#include <cstdint>
#include <ctime>
#include <string>
#include <vector>
#include <fstream>

/// Lecroy parser namespace
namespace lecroyparser {

/// Custom tm structure
struct time_stamp
{
    double  seconds;            ///< Seconds. [0-59]
    uint8_t minutes;            ///< Minutes. [0-59]
    uint8_t hours;              ///< Hours. [0-23]
    uint8_t days;               ///< Days. [1-31]
    uint8_t months;             ///< Months. [1-12]
    uint16_t year;              ///< Months. [0-16000]
    uint16_t unused;            ///< empty
} __attribute__((packed));///< packed struct

/// Template content of the LECROY_2_3 binary file. Adapted from https://github.com/sigrokproject/libsigrok/blob/master/src/hardware/lecroy-xstream/protocol.c
struct lecroy_wavedesc_2_3 {
    char descriptor_name[16];   ///< the first 8 chars will be WAVEDESC
    char template_name[16];     ///< LECROY_2_3
    uint16_t comm_type; ///< chosen by remote command COMM_FORMAT, 0: byte, 1: word
    uint16_t comm_order;///< 0: HIFIRST, 1: LOWFIRST, i.e. 1 - little endian
    uint32_t wave_descriptor_length;///< length in bytes of block WAVEDESC
    uint32_t user_text_len;///< length in bytes of block USERTEXT
    uint32_t res_desc1;    ///< an expansion entry is reserved
    uint32_t trigtime_array_length; ///< length in bytes of TRIGTIME array
    uint32_t ris_time1_array_length;///< length in bytes of RIS_TIME array
    uint32_t res_array1;            ///< an expansion entry is reserved
    uint32_t wave_array1_length;///< length in bytes of 1st simple data array. In transmitted waveform, represent the number of transmitted bytes in accordance with the NP parameter of the WFSU remote command and the used format (see COMM_TYPE).
    uint32_t wave_array2_length;///< length in bytes of 2nd simple data array
    uint32_t res_array2_length; ///< an expansion entry is reserved
    uint32_t res_array3_length; ///< an expansion entry is reserved
    char instrument_name[16];   ///< instrument name
    uint32_t instrument_number; ///< instrument number
    char trace_label[16];       ///< identifies the waveform
    uint16_t reserved1;         ///< expansion entry 1
    uint16_t reserved2;         ///< expansion entry 2
    uint32_t wave_array_count;  ///< number of data points in the data array. If there are two data arrays (FFT or Extrema), this number applies to each array separately.
    uint32_t points_per_screen; ///< nominal number of data points on the screen
    uint32_t first_valid_point; ///< count of number of points to skip before first good point. Equals 0 for normal waveforms.
    uint32_t last_valid_point;  ///< index of last good data point in record before padding (blanking) was started. LAST_VALID_POINT = WAVE_ARRAY_COUNT-1 except for aborted sequence and rollmode acquisitions
    uint32_t first_point;       ///< for input and output, indicates the offset relative to the beginning of the trace buffer. Value is the same as the FP parameter of the WFSU remote command.
    uint32_t sparsing_factor;   ///< for input and output, indicates the sparsing into the transmitted data block. Value is the same as the SP parameter of the WFSU remote command.
    uint32_t segment_index;     ///< for input and output, indicates the index of the transmitted segment. Value is the same as the SN parameter of the WFSU remote command.
    uint32_t subarray_count;    ///< for Sequence, acquired segment count, between 0 and NOM_SUBARRAY_COUNT
    uint32_t sweeps_per_acq;    ///< for Average or Extrema, number of sweeps accumulated, else 1
    uint16_t points_per_pair;   ///< for Peak Dectect waveforms (which always include data points in DATA_ARRAY_1 and min/max pairs in DATA_ARRAY_2). Value is the number of data points for each min/max pair.
    uint16_t pair_offset;       ///< for Peak Dectect waveforms only Value is the number of data points by which the first min/max pair in DATA_ARRAY_2 is offset relative to the first data value in DATA_ARRAY_1.
    float vertical_gain;        ///< gain
    float vertical_offset;      ///< to get floating values from raw data: VERTICAL_GAIN * data - VERTICAL_OFFSET
    float max_value;            ///< maximum allowed value. It corresponds to the upper edge of the grid.
    float min_value;            ///< minimum allowed value. It corresponds to the lower edge of the grid.
    uint16_t nominal_bits;      ///< a measure of the intrinsic precision of the observation: ADC data is 8 bit averaged data is 10-12 bit, etc.
    uint16_t nom_subarray_count;///< for Sequence, nominal segment count, else 1
    float horiz_interval;       ///< sampling interval for time domain waveforms
    double horiz_offset;        ///< trigger offset for the first sweep of the trigger, seconds between the trigger and the first data point
    double pixel_offset;        ///< needed to know how to display the waveform
    char vertunit[48];          ///< units of the vertical axis
    char horunit[48];           ///< units of the horizontal axis
    float horiz_uncertainty;    ///< uncertainty from one acquisition to the next, of the horizontal offset in seconds
    time_stamp trigger_time;    ///< time of the trigger
    float acq_duration;         ///< duration of the acquisition (in sec) in multi-trigger waveforms (e.g. sequence, RIS, or averaging)
    uint16_t record_type;       ///< RECORD_TYPE: enum
    uint16_t processing_done;   ///< PROCESSING_DONE: enum
    uint16_t reserved5;         ///< expansion entry
    uint16_t ris_sweeps;        ///< for RIS, the number of sweeps else 1
    uint16_t timebase;          ///< TIMEBASE: enum
    uint16_t vert_coupling;     ///< VERT_COUPLING: enum
    float probe_att;            ///< float
    uint16_t fixed_vert_gain;   ///< FIXED_VERT_GAIN: enum
    uint16_t bandwidth_limit;   ///< BANDWIDTH_LIMIT: enum
    float vertical_vernier;     ///< float
    float acq_vert_offset;      ///< float
    uint16_t wave_source;       ///< WAVE_SOURCE: enum
} __attribute__((packed));///< packed struct

/*
/// scope coupling options
static const char *scope_coupling_options[] = {
    "A1M", ///< AC with 1 MOhm termination
    "D50", ///< DC with 50 Ohm termination
    "D1M", ///< DC with 1 MOhm termination
    "GND",
    "OVL",
};

/// scope trigger slopes
static const char *scope_trigger_slopes[] = {
    "POS", "NEG",
};

/// scope trigger sources
static const char *scope_trigger_sources[] = {
    "C1", "C2", "C3", "C4", "LINE", "EXT",
};

/// scope timebases
static const uint64_t scope_timebases[][2] = {
    // picoseconds
    { 20, UINT64_C(1000000000000) },
    { 50, UINT64_C(1000000000000) },
    { 100, UINT64_C(1000000000000) },
    { 200, UINT64_C(1000000000000) },
    { 500, UINT64_C(1000000000000) },
    // nanoseconds
    { 1, 1000000000 },
    { 2, 1000000000 },
    { 5, 1000000000 },
    { 10, 1000000000 },
    { 20, 1000000000 },
    { 50, 1000000000 },
    { 100, 1000000000 },
    { 200, 1000000000 },
    { 500, 1000000000 },
    // microseconds
    { 1, 1000000 },
    { 2, 1000000 },
    { 5, 1000000 },
    { 10, 1000000 },
    { 20, 1000000 },
    { 50, 1000000 },
    { 100, 1000000 },
    { 200, 1000000 },
    { 500, 1000000 },
    // milliseconds
    { 1, 1000 },
    { 2, 1000 },
    { 5, 1000 },
    { 10, 1000 },
    { 20, 1000 },
    { 50, 1000 },
    { 100, 1000 },
    { 200, 1000 },
    { 500, 1000 },
    // seconds
    { 1, 1 },
    { 2, 1 },
    { 5, 1 },
    { 10, 1 },
    { 20, 1 },
    { 50, 1 },
    { 100, 1 },
    { 200, 1 },
    { 500, 1 },
    { 1000, 1 },
};

/// voltage division scales
static const uint64_t scope_vdivs[][2] = {
    // millivolts
    { 1, 1000 },
    { 2, 1000 },
    { 5, 1000 },
    { 10, 1000 },
    { 20, 1000 },
    { 50, 1000 },
    { 100, 1000 },
    { 200, 1000 },
    { 500, 1000 },
    // volts
    { 1, 1 },
    { 2, 1 },
    { 5, 1 },
    { 10, 1 },
    { 20, 1 },
    { 50, 1 },
};

/// channel names
static const char *scope_analog_channel_names[] = {
    "CH1", "CH2", "CH3", "CH4",
};
*/
/**
 * @brief read a trc file and store the header and waveform in the arguments passed as reference
 * @param fileName the name of the file to open
 * @param header where the header will be stored
 * @param waveform the pointer to the array where the waveform will be stored. Memory is allocated automatically if nullptr, and then object is owned by the caller, responsible of its deletion. If waveform is not nullptr, then it is assumed that user already preallocated the correct amount of memory, no checks are performed.
 * @return false in case an error is found, true otherwise
 * @return
 */
static bool read(const std::string fileName, lecroy_wavedesc_2_3& header, int16_t*& waveform) {
//static bool read(const std::string fileName, lecroy_wavedesc_2_3& header, long int*& waveform) {
    std::ifstream f(fileName, std::ios::binary);
    if(!f.good())
    {
        return false;
    }
    f.seekg (11, std::ios::beg);// skip magic bytes
    f.read((char*)&header, sizeof(header));
    if(header.wave_descriptor_length!=sizeof(header))
    {
        return header.wave_descriptor_length;
    }
    f.seekg(header.user_text_len+header.trigtime_array_length+header.ris_time1_array_length, std::ios::cur);
	if(!waveform) waveform = new int16_t[header.wave_array_count];
	//if(!waveform) waveform = new long int[header.wave_array_count];
    f.read((char*)waveform, header.wave_array1_length);
    return true;
}

}
