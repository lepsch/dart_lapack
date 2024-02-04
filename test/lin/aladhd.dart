      void aladhd(IOUNIT, PATH ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                IOUNIT;
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               CORZ, SORD;
      String             C1, C3;
      String             P2;
      String             SYM;
      // ..
      // .. External Functions ..
      //- bool               lsame, LSAMEN;
      // EXTERNAL lsame, LSAMEN
      // ..
      // .. Executable Statements ..

      if (IOUNIT <= 0) return;
      C1 = PATH( 1: 1 );
      C3 = PATH( 3: 3 );
      P2 = PATH( 2: 3 );
      SORD = lsame( C1, 'S' ) || lsame( C1, 'D' );
      CORZ = lsame( C1, 'C' ) || lsame( C1, 'Z' );
      if( !( SORD || CORZ ) ) return;

      if ( LSAMEN( 2, P2, 'GE' ) ) {

         // GE: General dense

         WRITE( IOUNIT, FMT = 9999 )PATH;
         WRITE( IOUNIT, FMT = '( '' Matrix types:'' )' );
         WRITE( IOUNIT, FMT = 9989 );
         WRITE( IOUNIT, FMT = '( '' Test ratios:'' )' );
         WRITE( IOUNIT, FMT = 9981 )1;
         WRITE( IOUNIT, FMT = 9980 )2;
         WRITE( IOUNIT, FMT = 9979 )3;
         WRITE( IOUNIT, FMT = 9978 )4;
         WRITE( IOUNIT, FMT = 9977 )5;
         WRITE( IOUNIT, FMT = 9976 )6;
         WRITE( IOUNIT, FMT = 9972 )7;
         WRITE( IOUNIT, FMT = '( '' Messages:'' )' );

      } else if ( LSAMEN( 2, P2, 'GB' ) ) {

         // GB: General band

         WRITE( IOUNIT, FMT = 9998 )PATH;
         WRITE( IOUNIT, FMT = '( '' Matrix types:'' )' );
         WRITE( IOUNIT, FMT = 9988 );
         WRITE( IOUNIT, FMT = '( '' Test ratios:'' )' );
         WRITE( IOUNIT, FMT = 9981 )1;
         WRITE( IOUNIT, FMT = 9980 )2;
         WRITE( IOUNIT, FMT = 9979 )3;
         WRITE( IOUNIT, FMT = 9978 )4;
         WRITE( IOUNIT, FMT = 9977 )5;
         WRITE( IOUNIT, FMT = 9976 )6;
         WRITE( IOUNIT, FMT = 9972 )7;
         WRITE( IOUNIT, FMT = '( '' Messages:'' )' );

      } else if ( LSAMEN( 2, P2, 'GT' ) ) {

         // GT: General tridiagonal

         WRITE( IOUNIT, FMT = 9997 )PATH;
         WRITE( IOUNIT, FMT = 9987 );
         WRITE( IOUNIT, FMT = '( '' Test ratios:'' )' );
         WRITE( IOUNIT, FMT = 9981 )1;
         WRITE( IOUNIT, FMT = 9980 )2;
         WRITE( IOUNIT, FMT = 9979 )3;
         WRITE( IOUNIT, FMT = 9978 )4;
         WRITE( IOUNIT, FMT = 9977 )5;
         WRITE( IOUNIT, FMT = 9976 )6;
         WRITE( IOUNIT, FMT = '( '' Messages:'' )' );

      } else if ( LSAMEN( 2, P2, 'PO' ) || LSAMEN( 2, P2, 'PP' ) || LSAMEN( 2, P2, 'PS' ) ) {

         // PO: Positive definite full
         // PS: Positive definite full
         // PP: Positive definite packed

         if ( SORD ) {
            SYM = 'Symmetric';
         } else {
            SYM = 'Hermitian';
         }
         if ( lsame( C3, 'O' ) ) {
            WRITE( IOUNIT, FMT = 9996 )PATH, SYM;
         } else {
            WRITE( IOUNIT, FMT = 9995 )PATH, SYM;
         }
         WRITE( IOUNIT, FMT = '( '' Matrix types:'' )' );
         WRITE( IOUNIT, FMT = 9985 )PATH;
         WRITE( IOUNIT, FMT = '( '' Test ratios:'' )' );
         WRITE( IOUNIT, FMT = 9975 )1;
         WRITE( IOUNIT, FMT = 9980 )2;
         WRITE( IOUNIT, FMT = 9979 )3;
         WRITE( IOUNIT, FMT = 9978 )4;
         WRITE( IOUNIT, FMT = 9977 )5;
         WRITE( IOUNIT, FMT = 9976 )6;
         WRITE( IOUNIT, FMT = '( '' Messages:'' )' );

      } else if ( LSAMEN( 2, P2, 'PB' ) ) {

         // PB: Positive definite band

         if ( SORD ) {
            WRITE( IOUNIT, FMT = 9994 )PATH, 'Symmetric';
         } else {
            WRITE( IOUNIT, FMT = 9994 )PATH, 'Hermitian';
         }
         WRITE( IOUNIT, FMT = '( '' Matrix types:'' )' );
         WRITE( IOUNIT, FMT = 9984 )PATH;
         WRITE( IOUNIT, FMT = '( '' Test ratios:'' )' );
         WRITE( IOUNIT, FMT = 9975 )1;
         WRITE( IOUNIT, FMT = 9980 )2;
         WRITE( IOUNIT, FMT = 9979 )3;
         WRITE( IOUNIT, FMT = 9978 )4;
         WRITE( IOUNIT, FMT = 9977 )5;
         WRITE( IOUNIT, FMT = 9976 )6;
         WRITE( IOUNIT, FMT = '( '' Messages:'' )' );

      } else if ( LSAMEN( 2, P2, 'PT' ) ) {

         // PT: Positive definite tridiagonal

         if ( SORD ) {
            WRITE( IOUNIT, FMT = 9993 )PATH, 'Symmetric';
         } else {
            WRITE( IOUNIT, FMT = 9993 )PATH, 'Hermitian';
         }
         WRITE( IOUNIT, FMT = 9986 );
         WRITE( IOUNIT, FMT = '( '' Test ratios:'' )' );
         WRITE( IOUNIT, FMT = 9973 )1;
         WRITE( IOUNIT, FMT = 9980 )2;
         WRITE( IOUNIT, FMT = 9979 )3;
         WRITE( IOUNIT, FMT = 9978 )4;
         WRITE( IOUNIT, FMT = 9977 )5;
         WRITE( IOUNIT, FMT = 9976 )6;
         WRITE( IOUNIT, FMT = '( '' Messages:'' )' );

      } else if ( LSAMEN( 2, P2, 'SY' ) || LSAMEN( 2, P2, 'SP' ) ) {

         // SY: Symmetric indefinite full
             // with partial (Bunch-Kaufman) pivoting algorithm
         // SP: Symmetric indefinite packed
             // with partial (Bunch-Kaufman) pivoting algorithm

         if ( lsame( C3, 'Y' ) ) {
            WRITE( IOUNIT, FMT = 9992 )PATH, 'Symmetric';
         } else {
            WRITE( IOUNIT, FMT = 9991 )PATH, 'Symmetric';
         }
         WRITE( IOUNIT, FMT = '( '' Matrix types:'' )' );
         if ( SORD ) {
            WRITE( IOUNIT, FMT = 9983 );
         } else {
            WRITE( IOUNIT, FMT = 9982 );
         }
         WRITE( IOUNIT, FMT = '( '' Test ratios:'' )' );
         WRITE( IOUNIT, FMT = 9974 )1;
         WRITE( IOUNIT, FMT = 9980 )2;
         WRITE( IOUNIT, FMT = 9979 )3;
         WRITE( IOUNIT, FMT = 9977 )4;
         WRITE( IOUNIT, FMT = 9978 )5;
         WRITE( IOUNIT, FMT = 9976 )6;
         WRITE( IOUNIT, FMT = '( '' Messages:'' )' );

      } else if ( LSAMEN( 2, P2, 'SR' ) || LSAMEN( 2, P2, 'SK') ) {

         // SR: Symmetric indefinite full,
             // with rook (bounded Bunch-Kaufman) pivoting algorithm

         // SK: Symmetric indefinite full,
             // with rook (bounded Bunch-Kaufman) pivoting algorithm,
             // ( new storage format for factors:
               // L and diagonal of D is stored in A,
               // subdiagonal of D is stored in E )

         WRITE( IOUNIT, FMT = 9992 )PATH, 'Symmetric';

         WRITE( IOUNIT, FMT = '( '' Matrix types:'' )' );
         if ( SORD ) {
            WRITE( IOUNIT, FMT = 9983 );
         } else {
            WRITE( IOUNIT, FMT = 9982 );
         }

         WRITE( IOUNIT, FMT = '( '' Test ratios:'' )' );
         WRITE( IOUNIT, FMT = 9974 )1;
         WRITE( IOUNIT, FMT = 9980 )2;
         WRITE( IOUNIT, FMT = 9979 )3;
         WRITE( IOUNIT, FMT = '( '' Messages:'' )' );

      } else if ( LSAMEN( 2, P2, 'HA' ) ) {

         // HA: Hermitian
             // Aasen algorithm
         WRITE( IOUNIT, FMT = 9971 )PATH, 'Hermitian';

         WRITE( IOUNIT, FMT = '( '' Matrix types:'' )' );
         WRITE( IOUNIT, FMT = 9983 );

         WRITE( IOUNIT, FMT = '( '' Test ratios:'' )' );
         WRITE( IOUNIT, FMT = 9974 )1;
         WRITE( IOUNIT, FMT = 9980 )2;
         WRITE( IOUNIT, FMT = 9979 )3;
         WRITE( IOUNIT, FMT = 9977 )4;
         WRITE( IOUNIT, FMT = 9978 )5;
         WRITE( IOUNIT, FMT = 9976 )6;
         WRITE( IOUNIT, FMT = '( '' Messages:'' )' );

       } else if ( LSAMEN( 2, P2, 'HE' ) || LSAMEN( 2, P2, 'HP' ) ) {

         // HE: Hermitian indefinite full
             // with partial (Bunch-Kaufman) pivoting algorithm
         // HP: Hermitian indefinite packed
             // with partial (Bunch-Kaufman) pivoting algorithm

         if ( lsame( C3, 'E' ) ) {
            WRITE( IOUNIT, FMT = 9992 )PATH, 'Hermitian';
         } else {
            WRITE( IOUNIT, FMT = 9991 )PATH, 'Hermitian';
         }

         WRITE( IOUNIT, FMT = '( '' Matrix types:'' )' );
         WRITE( IOUNIT, FMT = 9983 );

         WRITE( IOUNIT, FMT = '( '' Test ratios:'' )' );
         WRITE( IOUNIT, FMT = 9974 )1;
         WRITE( IOUNIT, FMT = 9980 )2;
         WRITE( IOUNIT, FMT = 9979 )3;
         WRITE( IOUNIT, FMT = 9977 )4;
         WRITE( IOUNIT, FMT = 9978 )5;
         WRITE( IOUNIT, FMT = 9976 )6;
         WRITE( IOUNIT, FMT = '( '' Messages:'' )' );

      } else if ( LSAMEN( 2, P2, 'HR' ) || LSAMEN( 2, P2, 'HK' ) ) {

         // HR: Hermitian indefinite full,
             // with rook (bounded Bunch-Kaufman) pivoting algorithm

         // HK: Hermitian indefinite full,
             // with rook (bounded Bunch-Kaufman) pivoting algorithm,
             // ( new storage format for factors:
               // L and diagonal of D is stored in A,
               // subdiagonal of D is stored in E )

         WRITE( IOUNIT, FMT = 9992 )PATH, 'Hermitian';

         WRITE( IOUNIT, FMT = '( '' Matrix types:'' )' );
         WRITE( IOUNIT, FMT = 9983 );

         WRITE( IOUNIT, FMT = '( '' Test ratios:'' )' );
         WRITE( IOUNIT, FMT = 9974 )1;
         WRITE( IOUNIT, FMT = 9980 )2;
         WRITE( IOUNIT, FMT = 9979 )3;
         WRITE( IOUNIT, FMT = '( '' Messages:'' )' );

      } else {

         // Print error message if no header is available.

         WRITE( IOUNIT, FMT = 9990 )PATH;
      }

      // First line of header

 9999 FORMAT( / 1X, A3, ' drivers:  General dense matrices' );
 9998 FORMAT( / 1X, A3, ' drivers:  General band matrices' );
 9997 FORMAT( / 1X, A3, ' drivers:  General tridiagonal' );
 9996 FORMAT( / 1X, A3, ' drivers:  ', A9, ' positive definite matrices' );
 9995 FORMAT( / 1X, A3, ' drivers:  ', A9, ' positive definite packed matrices' );
 9994 FORMAT( / 1X, A3, ' drivers:  ', A9, ' positive definite band matrices' );
 9993 FORMAT( / 1X, A3, ' drivers:  ', A9, ' positive definite tridiagonal' );
 9971 FORMAT( / 1X, A3, ' drivers:  ', A9, ' indefinite matrices', ', "Aasen" Algorithm' );
 9992 FORMAT( / 1X, A3, ' drivers:  ', A9, ' indefinite matrices', ', "rook" (bounded Bunch-Kaufman) pivoting' );
 9991 FORMAT( / 1X, A3, ' drivers:  ', A9, ' indefinite packed matrices', ', partial (Bunch-Kaufman) pivoting' );
 9891 FORMAT( / 1X, A3, ' drivers:  ', A9, ' indefinite packed matrices', ', "rook" (bounded Bunch-Kaufman) pivoting' );
 9990 FORMAT( / 1X, A3, ':  No header available' );

      // GE matrix types

 9989 FORMAT( 4X, '1. Diagonal', 24X, '7. Last n/2 columns zero', / 4X, '2. Upper triangular', 16X, '8. Random, CNDNUM = sqrt(0.1/EPS)', / 4X, '3. Lower triangular', 16X, '9. Random, CNDNUM = 0.1/EPS', / 4X, '4. Random, CNDNUM = 2', 13X, '10. Scaled near underflow', / 4X, '5. First column zero', 14X, '11. Scaled near overflow', / 4X, '6. Last column zero' );

      // GB matrix types

 9988 FORMAT( 4X, '1. Random, CNDNUM = 2', 14X, '5. Random, CNDNUM = sqrt(0.1/EPS)', / 4X, '2. First column zero', 15X, '6. Random, CNDNUM = 0.1/EPS', / 4X, '3. Last column zero', 16X, '7. Scaled near underflow', / 4X, '4. Last n/2 columns zero', 11X, '8. Scaled near overflow' );

      // GT matrix types

 9987 FORMAT( ' Matrix types (1-6 have specified condition numbers):', / 4X, '1. Diagonal', 24X, '7. Random, unspecified CNDNUM', / 4X, '2. Random, CNDNUM = 2', 14X, '8. First column zero', / 4X, '3. Random, CNDNUM = sqrt(0.1/EPS)', 2X, '9. Last column zero', / 4X, '4. Random, CNDNUM = 0.1/EPS', 7X, '10. Last n/2 columns zero', / 4X, '5. Scaled near underflow', 10X, '11. Scaled near underflow', / 4X, '6. Scaled near overflow', 11X, '12. Scaled near overflow' );

      // PT matrix types

 9986 FORMAT( ' Matrix types (1-6 have specified condition numbers):', / 4X, '1. Diagonal', 24X, '7. Random, unspecified CNDNUM', / 4X, '2. Random, CNDNUM = 2', 14X, '8. First row and column zero', / 4X, '3. Random, CNDNUM = sqrt(0.1/EPS)', 2X, '9. Last row and column zero', / 4X, '4. Random, CNDNUM = 0.1/EPS', 7X, '10. Middle row and column zero', / 4X, '5. Scaled near underflow', 10X, '11. Scaled near underflow', / 4X, '6. Scaled near overflow', 11X, '12. Scaled near overflow' );

      // PO, PP matrix types

 9985 FORMAT( 4X, '1. Diagonal', 24X, '6. Random, CNDNUM = sqrt(0.1/EPS)', / 4X, '2. Random, CNDNUM = 2', 14X, '7. Random, CNDNUM = 0.1/EPS', / 3X, '*3. First row and column zero', 7X, '8. Scaled near underflow', / 3X, '*4. Last row and column zero', 8X, '9. Scaled near overflow', / 3X, '*5. Middle row and column zero', / 3X, '(* - tests error exits from ', A3, 'TRF, no test ratios are computed)' );

      // PB matrix types

 9984 FORMAT( 4X, '1. Random, CNDNUM = 2', 14X, '5. Random, CNDNUM = sqrt(0.1/EPS)', / 3X, '*2. First row and column zero', 7X, '6. Random, CNDNUM = 0.1/EPS', / 3X, '*3. Last row and column zero', 8X, '7. Scaled near underflow', / 3X, '*4. Middle row and column zero', 6X, '8. Scaled near overflow', / 3X, '(* - tests error exits from ', A3, 'TRF, no test ratios are computed)' );

      // SSY, SSP, CHE, CHP matrix types

 9983 FORMAT( 4X, '1. Diagonal', 24X, '6. Last n/2 rows and columns zero', / 4X, '2. Random, CNDNUM = 2', 14X, '7. Random, CNDNUM = sqrt(0.1/EPS)', / 4X, '3. First row and column zero', 7X, '8. Random, CNDNUM = 0.1/EPS', / 4X, '4. Last row and column zero', 8X, '9. Scaled near underflow', / 4X, '5. Middle row and column zero', 5X, '10. Scaled near overflow' );

      // CSY, CSP matrix types

 9982 FORMAT( 4X, '1. Diagonal', 24X, '7. Random, CNDNUM = sqrt(0.1/EPS)', / 4X, '2. Random, CNDNUM = 2', 14X, '8. Random, CNDNUM = 0.1/EPS', / 4X, '3. First row and column zero', 7X, '9. Scaled near underflow', / 4X, '4. Last row and column zero', 7X, '10. Scaled near overflow', / 4X, '5. Middle row and column zero', 5X, '11. Block diagonal matrix', / 4X, '6. Last n/2 rows and columns zero' );

      // Test ratios

 9981 FORMAT( 3X, I2, ': norm( L * U - A )  / ( N * norm(A) * EPS )' );
 9980 FORMAT( 3X, I2, ': norm( B - A * X )  / ', '( norm(A) * norm(X) * EPS )' );
 9979 FORMAT( 3X, I2, ': norm( X - XACT )   / ', '( norm(XACT) * CNDNUM * EPS )' );
 9978 FORMAT( 3X, I2, ': norm( X - XACT )   / ', '( norm(XACT) * (error bound) )' );
 9977 FORMAT( 3X, I2, ': (backward error)   / EPS' );
 9976 FORMAT( 3X, I2, ': RCOND * CNDNUM - 1.0' );
 9975 FORMAT( 3X, I2, ': norm( U'' * U - A ) / ( N * norm(A) * EPS )', ', or', / 7X, 'norm( L * L'' - A ) / ( N * norm(A) * EPS )' );
 9974 FORMAT( 3X, I2, ': norm( U*D*U'' - A ) / ( N * norm(A) * EPS )', ', or', / 7X, 'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )' );
 9973 FORMAT( 3X, I2, ': norm( U''*D*U - A ) / ( N * norm(A) * EPS )', ', or', / 7X, 'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )' );
 9972 FORMAT( 3X, I2, ': abs( WORK(1) - RPVGRW ) /', ' ( max( WORK(1), RPVGRW ) * EPS )' );

      return;
      }