      void aladhd(final int IOUNIT, final int PATH) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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

      if (IOUNIT <= 0) return;
      C1 = PATH( 1: 1 );
      C3 = PATH( 3: 3 );
      P2 = PATH( 2: 3 );
      SORD = lsame( C1, 'S' ) || lsame( C1, 'D' );
      CORZ = lsame( C1, 'C' ) || lsame( C1, 'Z' );
      if( !( SORD || CORZ ) ) return;

      if ( lsamen( 2, P2, 'GE' ) ) {

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

      } else if ( lsamen( 2, P2, 'GB' ) ) {

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

      } else if ( lsamen( 2, P2, 'GT' ) ) {

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

      } else if ( lsamen( 2, P2, 'PO' ) || lsamen( 2, P2, 'PP' ) || lsamen( 2, P2, 'PS' ) ) {

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

      } else if ( lsamen( 2, P2, 'PB' ) ) {

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

      } else if ( lsamen( 2, P2, 'PT' ) ) {

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

      } else if ( lsamen( 2, P2, 'SY' ) || lsamen( 2, P2, 'SP' ) ) {

         // SY: Symmetric indefinite full
         //     with partial (Bunch-Kaufman) pivoting algorithm
         // SP: Symmetric indefinite packed
         //     with partial (Bunch-Kaufman) pivoting algorithm

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

      } else if ( lsamen( 2, P2, 'SR' ) || lsamen( 2, P2, 'SK') ) {

         // SR: Symmetric indefinite full,
         //     with rook (bounded Bunch-Kaufman) pivoting algorithm

         // SK: Symmetric indefinite full,
         //     with rook (bounded Bunch-Kaufman) pivoting algorithm,
         //     ( new storage format for factors:
         //       L and diagonal of D is stored in A,
         //       subdiagonal of D is stored in E )

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

      } else if ( lsamen( 2, P2, 'HA' ) ) {

         // HA: Hermitian
         //     Aasen algorithm
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

       } else if ( lsamen( 2, P2, 'HE' ) || lsamen( 2, P2, 'HP' ) ) {

         // HE: Hermitian indefinite full
         //     with partial (Bunch-Kaufman) pivoting algorithm
         // HP: Hermitian indefinite packed
         //     with partial (Bunch-Kaufman) pivoting algorithm

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

      } else if ( lsamen( 2, P2, 'HR' ) || lsamen( 2, P2, 'HK' ) ) {

         // HR: Hermitian indefinite full,
         //     with rook (bounded Bunch-Kaufman) pivoting algorithm

         // HK: Hermitian indefinite full,
         //     with rook (bounded Bunch-Kaufman) pivoting algorithm,
         //     ( new storage format for factors:
         //       L and diagonal of D is stored in A,
         //       subdiagonal of D is stored in E )

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

 9999 FORMAT('\n ${.a3} drivers:  General dense matrices' );
 9998 FORMAT('\n ${.a3} drivers:  General band matrices' );
 9997 FORMAT('\n ${.a3} drivers:  General tridiagonal' );
 9996 FORMAT('\n ${.a3} drivers:  ${.a9} positive definite matrices' );
 9995 FORMAT('\n ${.a3} drivers:  ${.a9} positive definite packed matrices' );
 9994 FORMAT('\n ${.a3} drivers:  ${.a9} positive definite band matrices' );
 9993 FORMAT('\n ${.a3} drivers:  ${.a9} positive definite tridiagonal' );
 9971 FORMAT('\n ${.a3} drivers:  ${.a9} indefinite matrices, "Aasen" Algorithm' );
 9992 FORMAT('\n ${.a3} drivers:  ${.a9} indefinite matrices, "rook" (bounded Bunch-Kaufman) pivoting' );
 9991 FORMAT('\n ${.a3} drivers:  ${.a9} indefinite packed matrices, partial (Bunch-Kaufman) pivoting' );
 9891 FORMAT('\n ${.a3} drivers:  ${.a9} indefinite packed matrices, "rook" (bounded Bunch-Kaufman) pivoting' );
 9990 FORMAT('\n ${.a3}:  No header available' );

      // GE matrix types

 9989 FORMAT('    1. Diagonal${' ' * 24}7. Last n/2 columns zero\n${' ' * 4}2. Upper triangular${' ' * 16}8. Random, CNDNUM = sqrt(0.1/EPS)\n${' ' * 4}3. Lower triangular${' ' * 16}9. Random, CNDNUM = 0.1/EPS\n${' ' * 4}4. Random, CNDNUM = 2${' ' * 13}10. Scaled near underflow\n${' ' * 4}5. First column zero${' ' * 14}11. Scaled near overflow\n${' ' * 4}6. Last column zero' );

      // GB matrix types

 9988 FORMAT('    1. Random, CNDNUM = 2${' ' * 14}5. Random, CNDNUM = sqrt(0.1/EPS)\n${' ' * 4}2. First column zero${' ' * 15}6. Random, CNDNUM = 0.1/EPS\n${' ' * 4}3. Last column zero${' ' * 16}7. Scaled near underflow\n${' ' * 4}4. Last n/2 columns zero${' ' * 11}8. Scaled near overflow' );

      // GT matrix types

 9987 FORMAT( ' Matrix types (1-6 have specified condition numbers):\n${' ' * 4}1. Diagonal${' ' * 24}7. Random, unspecified CNDNUM\n${' ' * 4}2. Random, CNDNUM = 2${' ' * 14}8. First column zero\n${' ' * 4}3. Random, CNDNUM = sqrt(0.1/EPS)${' ' * 2}9. Last column zero\n${' ' * 4}4. Random, CNDNUM = 0.1/EPS${' ' * 7}10. Last n/2 columns zero\n${' ' * 4}5. Scaled near underflow${' ' * 10}11. Scaled near underflow\n${' ' * 4}6. Scaled near overflow${' ' * 11}12. Scaled near overflow' );

      // PT matrix types

 9986 FORMAT( ' Matrix types (1-6 have specified condition numbers):\n${' ' * 4}1. Diagonal${' ' * 24}7. Random, unspecified CNDNUM\n${' ' * 4}2. Random, CNDNUM = 2${' ' * 14}8. First row and column zero\n${' ' * 4}3. Random, CNDNUM = sqrt(0.1/EPS)${' ' * 2}9. Last row and column zero\n${' ' * 4}4. Random, CNDNUM = 0.1/EPS${' ' * 7}10. Middle row and column zero\n${' ' * 4}5. Scaled near underflow${' ' * 10}11. Scaled near underflow\n${' ' * 4}6. Scaled near overflow${' ' * 11}12. Scaled near overflow' );

      // PO, PP matrix types

 9985 FORMAT('    1. Diagonal${' ' * 24}6. Random, CNDNUM = sqrt(0.1/EPS)\n${' ' * 4}2. Random, CNDNUM = 2${' ' * 14}7. Random, CNDNUM = 0.1/EPS\n${' ' * 3}*3. First row and column zero${' ' * 7}8. Scaled near underflow\n${' ' * 3}*4. Last row and column zero${' ' * 8}9. Scaled near overflow\n${' ' * 3}*5. Middle row and column zero\n${' ' * 3}(* - tests error exits from ${.a3}TRF, no test ratios are computed)' );

      // PB matrix types

 9984 FORMAT('    1. Random, CNDNUM = 2${' ' * 14}5. Random, CNDNUM = sqrt(0.1/EPS)\n${' ' * 3}*2. First row and column zero${' ' * 7}6. Random, CNDNUM = 0.1/EPS\n${' ' * 3}*3. Last row and column zero${' ' * 8}7. Scaled near underflow\n${' ' * 3}*4. Middle row and column zero${' ' * 6}8. Scaled near overflow\n${' ' * 3}(* - tests error exits from ${.a3}TRF, no test ratios are computed)' );

      // SSY, SSP, CHE, CHP matrix types

 9983 FORMAT('    1. Diagonal${' ' * 24}6. Last n/2 rows and columns zero\n${' ' * 4}2. Random, CNDNUM = 2${' ' * 14}7. Random, CNDNUM = sqrt(0.1/EPS)\n${' ' * 4}3. First row and column zero${' ' * 7}8. Random, CNDNUM = 0.1/EPS\n${' ' * 4}4. Last row and column zero${' ' * 8}9. Scaled near underflow\n${' ' * 4}5. Middle row and column zero${' ' * 5}10. Scaled near overflow' );

      // CSY, CSP matrix types

 9982 FORMAT('    1. Diagonal${' ' * 24}7. Random, CNDNUM = sqrt(0.1/EPS)\n${' ' * 4}2. Random, CNDNUM = 2${' ' * 14}8. Random, CNDNUM = 0.1/EPS\n${' ' * 4}3. First row and column zero${' ' * 7}9. Scaled near underflow\n${' ' * 4}4. Last row and column zero${' ' * 7}10. Scaled near overflow\n${' ' * 4}5. Middle row and column zero${' ' * 5}11. Block diagonal matrix\n${' ' * 4}6. Last n/2 rows and columns zero' );

      // Test ratios

 9981 FORMAT('${' ' * 3}${.i2}: norm( L * U - A )  / ( N * norm(A) * EPS )' );
 9980 FORMAT('${' ' * 3}${.i2}: norm( B - A * X )  / ( norm(A) * norm(X) * EPS )' );
 9979 FORMAT('${' ' * 3}${.i2}: norm( X - XACT )   / ( norm(XACT) * CNDNUM * EPS )' );
 9978 FORMAT('${' ' * 3}${.i2}: norm( X - XACT )   / ( norm(XACT) * (error bound) )' );
 9977 FORMAT('${' ' * 3}${.i2}: (backward error)   / EPS' );
 9976 FORMAT('${' ' * 3}${.i2}: RCOND * CNDNUM - 1.0' );
 9975 FORMAT('${' ' * 3}${.i2}: norm( U'' * U - A ) / ( N * norm(A) * EPS ), or\n${' ' * 7}norm( L * L'' - A ) / ( N * norm(A) * EPS )' );
 9974 FORMAT('${' ' * 3}${.i2}: norm( U*D*U'' - A ) / ( N * norm(A) * EPS ), or\n${' ' * 7}norm( L*D*L'' - A ) / ( N * norm(A) * EPS )' );
 9973 FORMAT('${' ' * 3}${.i2}: norm( U''*D*U - A ) / ( N * norm(A) * EPS ), or\n${' ' * 7}norm( L*D*L'' - A ) / ( N * norm(A) * EPS )' );
 9972 FORMAT('${' ' * 3}${.i2}: abs( WORK(1) - RPVGRW ) / ( max( WORK(1), RPVGRW ) * EPS )' );

      }
