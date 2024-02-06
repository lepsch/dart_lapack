      void main() {
// -- LAPACK test routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      int                IEEEOK;

      WRITE( 6, FMT = * ) 'We are about to check whether infinity arithmetic';
      WRITE( 6, FMT = * )'can be trusted.  If this test hangs, set';
      WRITE( 6, FMT = * ) 'ILAENV = 0 for ISPEC = 11 in LAPACK/SRC/ilaenv.f';

      IEEEOK = ilaenv( 11, 'ILAENV', 'N', 1, 2, 3, 4 );
      WRITE( 6, FMT = * );

      if ( IEEEOK == 0 ) {
         WRITE( 6, FMT = * ) 'Infinity arithmetic did not perform per the ieee spec';
      } else {
         WRITE( 6, FMT = * ) 'Infinity arithmetic performed as per the ieee spec.'          WRITE( 6, FMT = * ) 'However, this is not an exhaustive test and does not'          WRITE( 6, FMT = * ) 'guarantee that infinity arithmetic meets the', ' ieee spec.';
      }

      WRITE( 6, FMT = * );
      // ilaenv( 10, ...) checks both infinity and NaN arithmetic
      // infinity has already been checked so checking NaN now
      WRITE( 6, FMT = * ) 'We are about to check whether NaN arithmetic';
      WRITE( 6, FMT = * )'can be trusted.  If this test hangs, set';
      WRITE( 6, FMT = * ) 'ILAENV = 0 for ISPEC = 10 in LAPACK/SRC/ilaenv.f';
      IEEEOK = ilaenv( 10, 'ILAENV', 'N', 1, 2, 3, 4 );

      WRITE( 6, FMT = * );
      if ( IEEEOK == 0 ) {
         WRITE( 6, FMT = * ) 'NaN arithmetic did not perform per the ieee spec';
      } else {
         WRITE( 6, FMT = * )'NaN arithmetic performed as per the ieee', ' spec.'          WRITE( 6, FMT = * ) 'However, this is not an exhaustive test and does not'          WRITE( 6, FMT = * )'guarantee that NaN arithmetic meets the', ' ieee spec.';
      }
      WRITE( 6, FMT = * );

      }
