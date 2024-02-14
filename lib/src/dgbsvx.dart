import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgbsvx(final int FACT, final int TRANS, final int N, final int KL, final int KU, final int NRHS, final Matrix<double> AB_, final int LDAB, final Matrix<double> AFB_, final int LDAFB, final Array<int> IPIV_, final int EQUED, final int R, final int C, final Matrix<double> B_, final int LDB, final Matrix<double> X_, final int LDX, final int RCOND, final int FERR, final int BERR, final Array<double> _WORK_, final Array<int> IWORK_, final Box<int> INFO,) {
  final AB = AB_.dim();
  final AFB = AFB_.dim();
  final IPIV = IPIV_.dim();
  final B = B_.dim();
  final X = X_.dim();
  final _WORK = _WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             EQUED, FACT, TRANS;
      int                INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS;
      double             RCOND;
      int                IPIV( * ), IWORK( * );
      double             AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), BERR( * ), C( * ), FERR( * ), R( * ), WORK( * ), X( LDX, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
      String             NORM;
      int                I, INFEQU, J, J1, J2;
      double             AMAX, ANORM, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, RPVGRW, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLANGB, DLANTB;
      // EXTERNAL lsame, DLAMCH, DLANGB, DLANTB
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGBCON, DGBEQU, DGBRFS, DGBTRF, DGBTRS, DLACPY, DLAQGB, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      EQUIL = lsame( FACT, 'E' );
      NOTRAN = lsame( TRANS, 'N' );
      if ( NOFACT || EQUIL ) {
         EQUED = 'N';
         ROWEQU = false;
         COLEQU = false;
      } else {
         ROWEQU = lsame( EQUED, 'R' ) || lsame( EQUED, 'B' );
         COLEQU = lsame( EQUED, 'C' ) || lsame( EQUED, 'B' );
         SMLNUM = dlamch( 'Safe minimum' );
         BIGNUM = ONE / SMLNUM;
      }

      // Test the input parameters.

      if ( !NOFACT && !EQUIL && !lsame( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( KL < 0 ) {
         INFO = -4;
      } else if ( KU < 0 ) {
         INFO = -5;
      } else if ( NRHS < 0 ) {
         INFO = -6;
      } else if ( LDAB < KL+KU+1 ) {
         INFO = -8;
      } else if ( LDAFB < 2*KL+KU+1 ) {
         INFO = -10;
      } else if ( lsame( FACT, 'F' ) && !( ROWEQU || COLEQU || lsame( EQUED, 'N' ) ) ) {
         INFO = -12;
      } else {
         if ( ROWEQU ) {
            RCMIN = BIGNUM;
            RCMAX = ZERO;
            for (J = 1; J <= N; J++) { // 10
               RCMIN = min( RCMIN, R( J ) );
               RCMAX = max( RCMAX, R( J ) );
            } // 10
            if ( RCMIN <= ZERO ) {
               INFO = -13;
            } else if ( N > 0 ) {
               ROWCND = max( RCMIN, SMLNUM ) / min( RCMAX, BIGNUM );
            } else {
               ROWCND = ONE;
            }
         }
         if ( COLEQU && INFO == 0 ) {
            RCMIN = BIGNUM;
            RCMAX = ZERO;
            for (J = 1; J <= N; J++) { // 20
               RCMIN = min( RCMIN, C( J ) );
               RCMAX = max( RCMAX, C( J ) );
            } // 20
            if ( RCMIN <= ZERO ) {
               INFO = -14;
            } else if ( N > 0 ) {
               COLCND = max( RCMIN, SMLNUM ) / min( RCMAX, BIGNUM );
            } else {
               COLCND = ONE;
            }
         }
         if ( INFO == 0 ) {
            if ( LDB < max( 1, N ) ) {
               INFO = -16;
            } else if ( LDX < max( 1, N ) ) {
               INFO = -18;
            }
         }
      }

      if ( INFO != 0 ) {
         xerbla('DGBSVX', -INFO );
         return;
      }

      if ( EQUIL ) {

         // Compute row and column scalings to equilibrate the matrix A.

         dgbequ(N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFEQU );
         if ( INFEQU == 0 ) {

            // Equilibrate the matrix.

            dlaqgb(N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, EQUED );
            ROWEQU = lsame( EQUED, 'R' ) || lsame( EQUED, 'B' );
            COLEQU = lsame( EQUED, 'C' ) || lsame( EQUED, 'B' );
         }
      }

      // Scale the right hand side.

      if ( NOTRAN ) {
         if ( ROWEQU ) {
            for (J = 1; J <= NRHS; J++) { // 40
               for (I = 1; I <= N; I++) { // 30
                  B[I][J] = R( I )*B( I, J );
               } // 30
            } // 40
         }
      } else if ( COLEQU ) {
         for (J = 1; J <= NRHS; J++) { // 60
            for (I = 1; I <= N; I++) { // 50
               B[I][J] = C( I )*B( I, J );
            } // 50
         } // 60
      }

      if ( NOFACT || EQUIL ) {

         // Compute the LU factorization of the band matrix A.

         for (J = 1; J <= N; J++) { // 70
            J1 = max( J-KU, 1 );
            J2 = min( J+KL, N );
            dcopy(J2-J1+1, AB( KU+1-J+J1, J ), 1, AFB( KL+KU+1-J+J1, J ), 1 );
         } // 70

         dgbtrf(N, N, KL, KU, AFB, LDAFB, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {

            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            ANORM = ZERO;
            for (J = 1; J <= INFO; J++) { // 90
               for (I = max( KU+2-J, 1 ); I <= min( N+KU+1-J, KL+KU+1 ); I++) { // 80
                  ANORM = max( ANORM, ( AB( I, J ) ).abs() );
               } // 80
            } // 90
            RPVGRW = DLANTB( 'M', 'U', 'N', INFO, min( INFO-1, KL+KU ), AFB( max( 1, KL+KU+2-INFO ), 1 ), LDAFB, WORK );
            if ( RPVGRW == ZERO ) {
               RPVGRW = ONE;
            } else {
               RPVGRW = ANORM / RPVGRW;
            }
            WORK[1] = RPVGRW;
            RCOND = ZERO;
            return;
         }
      }

      // Compute the norm of the matrix A and the
      // reciprocal pivot growth factor RPVGRW.

      if ( NOTRAN ) {
         NORM = '1';
      } else {
         NORM = 'I';
      }
      ANORM = dlangb( NORM, N, KL, KU, AB, LDAB, WORK );
      RPVGRW = DLANTB( 'M', 'U', 'N', N, KL+KU, AFB, LDAFB, WORK );
      if ( RPVGRW == ZERO ) {
         RPVGRW = ONE;
      } else {
         RPVGRW = dlangb( 'M', N, KL, KU, AB, LDAB, WORK ) / RPVGRW;
      }

      // Compute the reciprocal of the condition number of A.

      dgbcon(NORM, N, KL, KU, AFB, LDAFB, IPIV, ANORM, RCOND, WORK, IWORK, INFO );

      // Compute the solution matrix X.

      dlacpy('Full', N, NRHS, B, LDB, X, LDX );
      dgbtrs(TRANS, N, KL, KU, NRHS, AFB, LDAFB, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      dgbrfs(TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO );

      // Transform the solution matrix X to a solution of the original
      // system.

      if ( NOTRAN ) {
         if ( COLEQU ) {
            for (J = 1; J <= NRHS; J++) { // 110
               for (I = 1; I <= N; I++) { // 100
                  X[I][J] = C( I )*X( I, J );
               } // 100
            } // 110
            for (J = 1; J <= NRHS; J++) { // 120
               FERR[J] = FERR( J ) / COLCND;
            } // 120
         }
      } else if ( ROWEQU ) {
         for (J = 1; J <= NRHS; J++) { // 140
            for (I = 1; I <= N; I++) { // 130
               X[I][J] = R( I )*X( I, J );
            } // 130
         } // 140
         for (J = 1; J <= NRHS; J++) { // 150
            FERR[J] = FERR( J ) / ROWCND;
         } // 150
      }

      // Set INFO = N+1 if the matrix is singular to working precision.

      if( RCOND < dlamch( 'Epsilon' ) ) INFO = N + 1;

      WORK[1] = RPVGRW;
      }
