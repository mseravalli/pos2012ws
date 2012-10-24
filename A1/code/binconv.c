#include <stdio.h>
#include <stdlib.h>

/**
 *
 */
int main ( int argc, char** argv ) {

    if ( argc < 3 ) {
        printf("Usage: ./binconv input_file output_file \n");
        exit( 1 );
    }
    
    char* in_path = argv[1]; 
    char* out_path = argv[2];

    FILE* in_file = fopen( in_path, "r" ); 
    FILE* out_file = fopen( out_path, "w" );
    
    if ( in_file == NULL ) {
        printf( "Error opening file: %s\n", in_path );
        exit( 1 );
    }

    // read dimensions

    int nintci; 
    int nintcf; 
    int nextci; 
    int nextcf; 

    fscanf( in_file, "%d", &nintci );
    fwrite( &nintci, sizeof( int ), 1, out_file );
    fscanf( in_file, "%d", &nintcf );
    fwrite( &nintcf, sizeof( int ), 1, out_file );
    fscanf( in_file, "%d", &nextci );
    fwrite( &nextci, sizeof( int ), 1, out_file );
    fscanf( in_file, "%d", &nextcf );
    fwrite( &nextcf, sizeof( int ), 1, out_file );

    int i;

    int temp_int;
    for ( i = nintci; i <= nintcf; i++ ) {
        fscanf( in_file, "%d", &temp_int );
        fwrite( &temp_int, sizeof( int ), 1, out_file );
        fscanf( in_file, "%d", &temp_int );
        fwrite( &temp_int, sizeof( int ), 1, out_file );
        fscanf( in_file, "%d", &temp_int );
        fwrite( &temp_int, sizeof( int ), 1, out_file );
        fscanf( in_file, "%d", &temp_int );
        fwrite( &temp_int, sizeof( int ), 1, out_file );
        fscanf( in_file, "%d", &temp_int );
        fwrite( &temp_int, sizeof( int ), 1, out_file );
        fscanf( in_file, "%d", &temp_int );
        fwrite( &temp_int, sizeof( int ), 1, out_file );
    }

    double temp_double;
    for ( i = nintci; i <= nintcf; i++ ) {
        fscanf( in_file, "%lf", &temp_double );
        fwrite( &temp_double, sizeof( double ), 1, out_file );
        fscanf( in_file, "%lf", &temp_double );
        fwrite( &temp_double, sizeof( double ), 1, out_file );
        fscanf( in_file, "%lf", &temp_double );
        fwrite( &temp_double, sizeof( double ), 1, out_file );
        fscanf( in_file, "%lf", &temp_double );
        fwrite( &temp_double, sizeof( double ), 1, out_file );
        fscanf( in_file, "%lf", &temp_double );
        fwrite( &temp_double, sizeof( double ), 1, out_file );
        fscanf( in_file, "%lf", &temp_double );
        fwrite( &temp_double, sizeof( double ), 1, out_file );
        fscanf( in_file, "%lf", &temp_double );
        fwrite( &temp_double, sizeof( double ), 1, out_file );
        fscanf( in_file, "%lf", &temp_double );
        fwrite( &temp_double, sizeof( double ), 1, out_file );
    }

    for ( i = nintci; i <= nintcf; i++ ) {
        fscanf( in_file, "%d", &temp_int );
        fwrite( &temp_int, sizeof( int ), 1, out_file );
    }

//  while ( !feof(in_file) ) {
//      fprintf( out_file, "%x ", next_int );
//      fflush( out_file );
//     
//      fscanf( in_file, "%d", &next_int );
//  }

    fclose(in_file);
    fclose(out_file);

    return 0;

}



