      SUBROUTINE ALAHDG( IOUNIT, PATH );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String            PATH;
      int               IOUNIT;
      // ..

*  =====================================================================

      // .. Local Scalars ..
      String            C2;
      int               ITYPE;
      // ..
      // .. External Functions ..
      bool              LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. Executable Statements ..

      if (IOUNIT <= 0) RETURN;
      C2 = PATH( 1: 3 );

      // First line describing matrices in this path

      if ( LSAMEN( 3, C2, 'GQR' ) ) {
         ITYPE = 1;
         WRITE( IOUNIT, FMT = 9991 )PATH;
      } else if ( LSAMEN( 3, C2, 'GRQ' ) ) {
         ITYPE = 2;
         WRITE( IOUNIT, FMT = 9992 )PATH;
      } else if ( LSAMEN( 3, C2, 'LSE' ) ) {
         ITYPE = 3;
         WRITE( IOUNIT, FMT = 9993 )PATH;
      } else if ( LSAMEN( 3, C2, 'GLM' ) ) {
         ITYPE = 4;
         WRITE( IOUNIT, FMT = 9994 )PATH;
      } else if ( LSAMEN( 3, C2, 'GSV' ) ) {
         ITYPE = 5;
         WRITE( IOUNIT, FMT = 9995 )PATH;
      } else if ( LSAMEN( 3, C2, 'CSD' ) ) {
         ITYPE = 6;
         WRITE( IOUNIT, FMT = 9996 )PATH;
      }

      // Matrix types

      WRITE( IOUNIT, FMT = 9999 )'Matrix types: ';

      if ( ITYPE == 1 ) {
         WRITE( IOUNIT, FMT = 9950 )1;
         WRITE( IOUNIT, FMT = 9952 )2;
         WRITE( IOUNIT, FMT = 9954 )3;
         WRITE( IOUNIT, FMT = 9955 )4;
         WRITE( IOUNIT, FMT = 9956 )5;
         WRITE( IOUNIT, FMT = 9957 )6;
         WRITE( IOUNIT, FMT = 9961 )7;
         WRITE( IOUNIT, FMT = 9962 )8;
      } else if ( ITYPE == 2 ) {
         WRITE( IOUNIT, FMT = 9951 )1;
         WRITE( IOUNIT, FMT = 9953 )2;
         WRITE( IOUNIT, FMT = 9954 )3;
         WRITE( IOUNIT, FMT = 9955 )4;
         WRITE( IOUNIT, FMT = 9956 )5;
         WRITE( IOUNIT, FMT = 9957 )6;
         WRITE( IOUNIT, FMT = 9961 )7;
         WRITE( IOUNIT, FMT = 9962 )8;
      } else if ( ITYPE == 3 ) {
         WRITE( IOUNIT, FMT = 9950 )1;
         WRITE( IOUNIT, FMT = 9952 )2;
         WRITE( IOUNIT, FMT = 9954 )3;
         WRITE( IOUNIT, FMT = 9955 )4;
         WRITE( IOUNIT, FMT = 9955 )5;
         WRITE( IOUNIT, FMT = 9955 )6;
         WRITE( IOUNIT, FMT = 9955 )7;
         WRITE( IOUNIT, FMT = 9955 )8;
      } else if ( ITYPE == 4 ) {
         WRITE( IOUNIT, FMT = 9951 )1;
         WRITE( IOUNIT, FMT = 9953 )2;
         WRITE( IOUNIT, FMT = 9954 )3;
         WRITE( IOUNIT, FMT = 9955 )4;
         WRITE( IOUNIT, FMT = 9955 )5;
         WRITE( IOUNIT, FMT = 9955 )6;
         WRITE( IOUNIT, FMT = 9955 )7;
         WRITE( IOUNIT, FMT = 9955 )8;
      } else if ( ITYPE == 5 ) {
         WRITE( IOUNIT, FMT = 9950 )1;
         WRITE( IOUNIT, FMT = 9952 )2;
         WRITE( IOUNIT, FMT = 9954 )3;
         WRITE( IOUNIT, FMT = 9955 )4;
         WRITE( IOUNIT, FMT = 9956 )5;
         WRITE( IOUNIT, FMT = 9957 )6;
         WRITE( IOUNIT, FMT = 9959 )7;
         WRITE( IOUNIT, FMT = 9960 )8;
      } else if ( ITYPE == 6 ) {
         WRITE( IOUNIT, FMT = 9963 )1;
         WRITE( IOUNIT, FMT = 9964 )2;
         WRITE( IOUNIT, FMT = 9965 )3;
      }

      // Tests performed

      WRITE( IOUNIT, FMT = 9999 )'Test ratios: ';

      if ( ITYPE == 1 ) {

         // GQR decomposition of rectangular matrices

         WRITE( IOUNIT, FMT = 9930 )1;
         WRITE( IOUNIT, FMT = 9931 )2;
         WRITE( IOUNIT, FMT = 9932 )3;
         WRITE( IOUNIT, FMT = 9933 )4;
      } else if ( ITYPE == 2 ) {

         // GRQ decomposition of rectangular matrices

         WRITE( IOUNIT, FMT = 9934 )1;
         WRITE( IOUNIT, FMT = 9935 )2;
         WRITE( IOUNIT, FMT = 9932 )3;
         WRITE( IOUNIT, FMT = 9933 )4;
      } else if ( ITYPE == 3 ) {

         // LSE Problem

         WRITE( IOUNIT, FMT = 9937 )1;
         WRITE( IOUNIT, FMT = 9938 )2;
      } else if ( ITYPE == 4 ) {

         // GLM Problem

         WRITE( IOUNIT, FMT = 9939 )1;
      } else if ( ITYPE == 5 ) {

         // GSVD

         WRITE( IOUNIT, FMT = 9940 )1;
         WRITE( IOUNIT, FMT = 9941 )2;
         WRITE( IOUNIT, FMT = 9942 )3;
         WRITE( IOUNIT, FMT = 9943 )4;
         WRITE( IOUNIT, FMT = 9944 )5;
      } else if ( ITYPE == 6 ) {

         // CSD

         WRITE( IOUNIT, FMT = 9910 );
         WRITE( IOUNIT, FMT = 9911 )1;
         WRITE( IOUNIT, FMT = 9912 )2;
         WRITE( IOUNIT, FMT = 9913 )3;
         WRITE( IOUNIT, FMT = 9914 )4;
         WRITE( IOUNIT, FMT = 9915 )5;
         WRITE( IOUNIT, FMT = 9916 )6;
         WRITE( IOUNIT, FMT = 9917 )7;
         WRITE( IOUNIT, FMT = 9918 )8;
         WRITE( IOUNIT, FMT = 9919 )9;
         WRITE( IOUNIT, FMT = 9920 );
         WRITE( IOUNIT, FMT = 9921 )10;
         WRITE( IOUNIT, FMT = 9922 )11;
         WRITE( IOUNIT, FMT = 9923 )12;
         WRITE( IOUNIT, FMT = 9924 )13;
         WRITE( IOUNIT, FMT = 9925 )14;
         WRITE( IOUNIT, FMT = 9926 )15;
      }

 9999 FORMAT( 1X, A );
 9991 FORMAT( / 1X, A3, ': GQR factorization of general matrices' );
 9992 FORMAT( / 1X, A3, ': GRQ factorization of general matrices' );
 9993 FORMAT( / 1X, A3, ': LSE Problem' );
 9994 FORMAT( / 1X, A3, ': GLM Problem' );
 9995 FORMAT( / 1X, A3, ': Generalized Singular Value Decomposition' );
 9996 FORMAT( / 1X, A3, ': CS Decomposition' );

 9950 FORMAT( 3X, I2, ': A-diagonal matrix  B-upper triangular' );
 9951 FORMAT( 3X, I2, ': A-diagonal matrix  B-lower triangular' );
 9952 FORMAT( 3X, I2, ': A-upper triangular B-upper triangular' );
 9953 FORMAT( 3X, I2, ': A-lower triangular B-diagonal triangular' );
 9954 FORMAT( 3X, I2, ': A-lower triangular B-upper triangular' );

 9955 FORMAT( 3X, I2, ': Random matrices cond(A)=100, cond(B)=10,' );

 9956 FORMAT( 3X, I2, ': Random matrices cond(A)= sqrt( 0.1/EPS ) ', 'cond(B)= sqrt( 0.1/EPS )' );
 9957 FORMAT( 3X, I2, ': Random matrices cond(A)= 0.1/EPS ', 'cond(B)= 0.1/EPS' );
 9959 FORMAT( 3X, I2, ': Random matrices cond(A)= sqrt( 0.1/EPS ) ', 'cond(B)=  0.1/EPS ' );
 9960 FORMAT( 3X, I2, ': Random matrices cond(A)= 0.1/EPS ', 'cond(B)=  sqrt( 0.1/EPS )' );

 9961 FORMAT( 3X, I2, ': Matrix scaled near underflow limit' );
 9962 FORMAT( 3X, I2, ': Matrix scaled near overflow limit' );
 9963 FORMAT( 3X, I2, ': Random orthogonal matrix (Haar measure)' );
 9964 FORMAT( 3X, I2, ': Nearly orthogonal matrix with uniformly ', 'distributed angles atan2( S, C ) in CS decomposition' );
 9965 FORMAT( 3X, I2, ': Random orthogonal matrix with clustered ', 'angles atan2( S, C ) in CS decomposition' );


      // GQR test ratio

 9930 FORMAT( 3X, I2, ': norm( R - Q'' * A ) / ( min( N, M )*norm( A )', '* EPS )' );
 9931 FORMAT( 3X, I2, ': norm( T * Z - Q'' * B )  / ( min(P,N)*norm(B)', '* EPS )' );
 9932 FORMAT( 3X, I2, ': norm( I - Q''*Q )   / ( N * EPS )' );
 9933 FORMAT( 3X, I2, ': norm( I - Z''*Z )   / ( P * EPS )' );

      // GRQ test ratio

 9934 FORMAT( 3X, I2, ': norm( R - A * Q'' ) / ( min( N,M )*norm(A) * ', 'EPS )' );
 9935 FORMAT( 3X, I2, ': norm( T * Q - Z'' * B )  / ( min( P,N ) * nor', 'm(B)*EPS )' );

      // LSE test ratio

 9937 FORMAT( 3X, I2, ': norm( A*x - c )  / ( norm(A)*norm(x) * EPS )' );
 9938 FORMAT( 3X, I2, ': norm( B*x - d )  / ( norm(B)*norm(x) * EPS )' );

      // GLM test ratio

 9939 FORMAT( 3X, I2, ': norm( d - A*x - B*y ) / ( (norm(A)+norm(B) )*', '(norm(x)+norm(y))*EPS )' );

      // GSVD test ratio

 9940 FORMAT( 3X, I2, ': norm( U'' * A * Q - D1 * R ) / ( min( M, N )*', 'norm( A ) * EPS )' );
 9941 FORMAT( 3X, I2, ': norm( V'' * B * Q - D2 * R ) / ( min( P, N )*', 'norm( B ) * EPS )' );
 9942 FORMAT( 3X, I2, ': norm( I - U''*U )   / ( M * EPS )' );
 9943 FORMAT( 3X, I2, ': norm( I - V''*V )   / ( P * EPS )' );
 9944 FORMAT( 3X, I2, ': norm( I - Q''*Q )   / ( N * EPS )' );

      // CSD test ratio

 9910 FORMAT( 3X, '2-by-2 CSD' );
 9911 FORMAT( 3X, I2, ': norm( U1'' * X11 * V1 - C ) / ( max(  P,  Q)', ' * max(norm(I-X''*X),EPS) )' );
 9912 FORMAT( 3X, I2, ': norm( U1'' * X12 * V2-(-S)) / ( max(  P,', 'M-Q) * max(norm(I-X''*X),EPS) )' );
 9913 FORMAT( 3X, I2, ': norm( U2'' * X21 * V1 - S ) / ( max(M-P,', '  Q) * max(norm(I-X''*X),EPS) )' );
 9914 FORMAT( 3X, I2, ': norm( U2'' * X22 * V2 - C ) / ( max(M-P,', 'M-Q) * max(norm(I-X''*X),EPS) )' );
 9915 FORMAT( 3X, I2, ': norm( I - U1''*U1 ) / (   P   * EPS )' );
 9916 FORMAT( 3X, I2, ': norm( I - U2''*U2 ) / ( (M-P) * EPS )' );
 9917 FORMAT( 3X, I2, ': norm( I - V1''*V1 ) / (   Q   * EPS )' );
 9918 FORMAT( 3X, I2, ': norm( I - V2''*V2 ) / ( (M-Q) * EPS )' );
 9919 FORMAT( 3X, I2, ': principal angle ordering ( 0 or ULP )' );
 9920 FORMAT( 3X, '2-by-1 CSD' );
 9921 FORMAT( 3X, I2, ': norm( U1'' * X11 * V1 - C ) / ( max(  P,  Q)', ' * max(norm(I-X''*X),EPS) )' );
 9922 FORMAT( 3X, I2, ': norm( U2'' * X21 * V1 - S ) / ( max(  M-P,', 'Q) * max(norm(I-X''*X),EPS) )' );
 9923 FORMAT( 3X, I2, ': norm( I - U1''*U1 ) / (   P   * EPS )' );
 9924 FORMAT( 3X, I2, ': norm( I - U2''*U2 ) / ( (M-P) * EPS )' );
 9925 FORMAT( 3X, I2, ': norm( I - V1''*V1 ) / (   Q   * EPS )' );
 9926 FORMAT( 3X, I2, ': principal angle ordering ( 0 or ULP )' );
      return;

      // End of ALAHDG

      }
