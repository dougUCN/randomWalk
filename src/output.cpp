#include <output.hpp>
#include <cmath>

void writeToHDF( std::string filename, std::vector<hdfOutputFormat> results , std::map<std::string, double> attributes)
{
    hsize_t nrecords = results.size();
    size_t dst_size =  sizeof( hdfOutputFormat );
    int *fill_data = NULL;
    hid_t      file_id;

    // Field information
    size_t dst_offset[NFIELDS] = {  HOFFSET( hdfOutputFormat,    particleNum),
                                    HOFFSET( hdfOutputFormat,    windowHits),
                                    HOFFSET( hdfOutputFormat,    totalSteps),
                                    HOFFSET( hdfOutputFormat,    location),
                                    HOFFSET( hdfOutputFormat,    status),
                                    HOFFSET( hdfOutputFormat,    cellRejections),
                                    HOFFSET( hdfOutputFormat,    cellExits),
                                    HOFFSET( hdfOutputFormat,    velocity),
                                    HOFFSET( hdfOutputFormat,    tstart),
                                    HOFFSET( hdfOutputFormat,    tend )};
    const char *field_names[NFIELDS]  = {   "particleNum",
                                            "windowHits",
                                            "totalSteps",
                                            "location",
                                            "status",
                                            "cellRejections",
                                            "cellExits",
                                            "velocity",
                                            "tstart",
                                            "tend"};

    // Initialize field types
    hid_t string_type = H5Tcopy( H5T_C_S1 );
    H5Tset_size( string_type, CHAR_COUNT );
    hid_t field_type[NFIELDS] = {    H5T_NATIVE_INT,
                            H5T_NATIVE_INT,
                            H5T_NATIVE_INT,
                            H5T_NATIVE_DOUBLE,
                            string_type,
                            H5T_NATIVE_INT,
                            H5T_NATIVE_INT,
                            H5T_NATIVE_DOUBLE,
                            H5T_NATIVE_DOUBLE,
                            H5T_NATIVE_DOUBLE};

    // Calculate file writing chunk size... this is kind of arbitrary since I couldn't
    // find a method to auto-calculate chunk size for you
    int chunk_size = std::ceil( nrecords / NFIELDS);
    double buffer;


    // Create file
    file_id = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    H5TBmake_table( "table", file_id, "table", NFIELDS, nrecords ,
                            dst_size, field_names, dst_offset, field_type,
                            chunk_size, fill_data, COMPRESSION, static_cast<void*>(results.data())  );

    // Write attributes
    for (auto const& it : attributes)
    {
        H5LTset_attribute_double(file_id, "/table", it.first.c_str(), &it.second, 1);
        H5LTget_attribute_double(file_id, "/table", it.first.c_str(), &buffer);
    }

    // Close everything else
    H5Tclose( string_type );
    H5Fclose( file_id );

}
