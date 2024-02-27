      import 'dart:math';

import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/ztrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zggrqf.dart';
import 'package:lapack/src/ztrtrs.dart';
import 'package:lapack/src/zunmqr.dart';
import 'package:lapack/src/zunmrq.dart';

void zgglse(final int M, final int N, final int P, final Matrix<Complex> A_,
      final int LDA, final Matrix<Complex> B_, final int LDB, final Array<Complex> C_,
      final Array<Complex> D_, final Array<Complex> X_, final Array<Complex> WORK_, final int LWORK,
      final Box<int> INFO,) {
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final WORK = WORK_.dim();
final C=C_.dim();
final D=D_.dim();
final X=X_.dim();
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               LQUERY;
      int                LOPT, LWKMIN, LWKOPT, MN, NB, NB1, NB2, NB3, NB4, NR;

      // Test the input parameters

      INFO.value = 0;
      MN = min( M, N );
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO.value = -1;
      } else if ( N < 0 ) {
         INFO.value = -2;
      } else if ( P < 0 || P > N || P < N-M ) {
         INFO.value = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO.value = -5;
      } else if ( LDB < max( 1, P ) ) {
         INFO.value = -7;
      }

      // Calculate workspace

      if ( INFO.value == 0) {
         if ( N == 0 ) {
            LWKMIN = 1;
            LWKOPT = 1;
         } else {
            NB1 = ilaenv( 1, 'ZGEQRF', ' ', M, N, -1, -1 );
            NB2 = ilaenv( 1, 'ZGERQF', ' ', M, N, -1, -1 );
            NB3 = ilaenv( 1, 'ZUNMQR', ' ', M, N, P, -1 );
            NB4 = ilaenv( 1, 'ZUNMRQ', ' ', M, N, P, -1 );
            NB = max( max(NB1, NB2), max(NB3, NB4) );
            LWKMIN = M + N + P;
            LWKOPT = P + MN + max( M, N )*NB;
         }
         WORK[1] = LWKOPT.toComplex();

         if ( LWORK < LWKMIN && !LQUERY ) {
            INFO.value = -12;
         }
      }

      if ( INFO.value != 0 ) {
         xerbla('ZGGLSE', -INFO.value );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Compute the GRQ factorization of matrices B and A:

      //        B*Q**H = (  0  T12 ) P   Z**H*A*Q**H = ( R11 R12 ) N-P
      //                    N-P  P                     (  0  R22 ) M+P-N
      //                                                  N-P  P

      // where T12 and R11 are upper triangular, and Q and Z are
      // unitary.

      zggrqf(P, M, N, B, LDB, WORK, A, LDA, WORK( P+1 ), WORK( P+MN+1 ), LWORK-P-MN, INFO );
      LOPT = WORK[ P+MN+1 ].toInt();

      // Update c = Z**H *c = ( c1 ) N-P
      //                   ( c2 ) M+P-N

      zunmqr('Left', 'Conjugate Transpose', M, 1, MN, A, LDA, WORK( P+1 ), C.asMatrix(max(1, M)), max( 1, M ), WORK( P+MN+1 ), LWORK-P-MN, INFO );
      LOPT = max( LOPT, WORK[ P+MN+1 ].toInt() );

      // Solve T12*x2 = d for x2

      if ( P > 0 ) {
         ztrtrs('Upper', 'No transpose', 'Non-unit', P, 1, B( 1, N-P+1 ), LDB, D.asMatrix(P), P, INFO );

         if ( INFO.value > 0 ) {
            INFO.value = 1;
            return;
         }

         // Put the solution in X

         zcopy(P, D, 1, X( N-P+1 ), 1 );

         // Update c1

         zgemv('No transpose', N-P, P, -Complex.one, A( 1, N-P+1 ), LDA, D, 1, Complex.one, C, 1 );
      }

      // Solve R11*x1 = c1 for x1

      if ( N > P ) {
         ztrtrs('Upper', 'No transpose', 'Non-unit', N-P, 1, A, LDA, C.asMatrix(N - P), N-P, INFO );

         if ( INFO.value > 0 ) {
            INFO.value = 2;
            return;
         }

         // Put the solutions in X

         zcopy(N-P, C, 1, X, 1 );
      }

      // Compute the residual vector:

      if ( M < N ) {
         NR = M + P - N;
         if (NR > 0) zgemv( 'No transpose', NR, N-M, -Complex.one, A( N-P+1, M+1 ), LDA, D( NR+1 ), 1, Complex.one, C( N-P+1 ), 1 );
      } else {
         NR = P;
      }
      if ( NR > 0 ) {
         ztrmv('Upper', 'No transpose', 'Non unit', NR, A( N-P+1, N-P+1 ), LDA, D, 1 );
         zaxpy(NR, -Complex.one, D, 1, C( N-P+1 ), 1 );
      }

      // Backward transformation x = Q**H*x

      zunmrq('Left', 'Conjugate Transpose', N, 1, P, B, LDB, WORK( 1 ), X.asMatrix(N), N, WORK( P+MN+1 ), LWORK-P-MN, INFO );
      WORK[1] = (P + MN + max( LOPT, WORK[ P+MN+1 ].toInt() )).toComplex();

      }
