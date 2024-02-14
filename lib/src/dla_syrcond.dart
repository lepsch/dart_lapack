import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      double dla_syrcond(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> AF_, final int LDAF, final Array<int> IPIV_, final int CMODE, final int C, final int INFO, final Array<double> _WORK_, final Array<int> IWORK_,) {
  final A = A_.dim();
  final AF = AF_.dim();
  final IPIV = IPIV_.dim();
  final _WORK = _WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                N, LDA, LDAF, INFO, CMODE;
      // ..
      // .. Array Arguments
      int                IWORK( * ), IPIV( * );
      double             A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      String             NORMIN;
      int                KASE, I, J;
      double             AINVNM, SMLNUM, TMP;
      bool               UP;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH;
      // EXTERNAL lsame, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, XERBLA, DSYTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX

      DLA_SYRCOND = 0.0;

      INFO = 0;
      if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( LDAF < max( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('DLA_SYRCOND', -INFO );
         return;
      }
      if ( N == 0 ) {
         DLA_SYRCOND = 1.0;
         return;
      }
      UP = false;
      if ( lsame( UPLO, 'U' ) ) UP = true;

      // Compute the equilibration matrix R such that
      // inv(R)*A*C has unit 1-norm.

      if ( UP ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            if ( CMODE == 1 ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( J, I ) * C( J ) );
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) * C( J ) );
               }
            } else if ( CMODE == 0 ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ( A( J, I ) ).abs();
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ( A( I, J ) ).abs();
               }
            } else {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( J, I ) / C( J ) );
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) / C( J ) );
               }
            }
            WORK[2*N+I] = TMP;
         }
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            if ( CMODE == 1 ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( I, J ) * C( J ) );
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) * C( J ) );
               }
            } else if ( CMODE == 0 ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ( A( I, J ) ).abs();
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ( A( J, I ) ).abs();
               }
            } else {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( I, J) / C( J ) );
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I) / C( J ) );
               }
            }
            WORK[2*N+I] = TMP;
         }
      }

      // Estimate the norm of inv(op(A)).

      SMLNUM = dlamch( 'Safe minimum' );
      AINVNM = 0.0;
      NORMIN = 'N';

      KASE = 0;
      // } // 10
      dlacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == 2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK[I] = WORK( I ) * WORK( 2*N+I );
            }

            if ( UP ) {
               dsytrs('U', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               dsytrs('L', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by inv(C).

            if ( CMODE == 1 ) {
               for (I = 1; I <= N; I++) {
                  WORK[I] = WORK( I ) / C( I );
               }
            } else if ( CMODE == -1 ) {
               for (I = 1; I <= N; I++) {
                  WORK[I] = WORK( I ) * C( I );
               }
            }
         } else {

            // Multiply by inv(C**T).

            if ( CMODE == 1 ) {
               for (I = 1; I <= N; I++) {
                  WORK[I] = WORK( I ) / C( I );
               }
            } else if ( CMODE == -1 ) {
               for (I = 1; I <= N; I++) {
                  WORK[I] = WORK( I ) * C( I );
               }
            }

            if ( UP ) {
               dsytrs('U', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               dsytrs('L', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK[I] = WORK( I ) * WORK( 2*N+I );
            }
         }

         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != 0.0) DLA_SYRCOND = ( 1.0 / AINVNM );

      }
