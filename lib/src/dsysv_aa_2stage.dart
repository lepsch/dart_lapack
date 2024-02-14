import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dsysv_aa_2stage(final int UPLO, final int N, final int NRHS, final Matrix<double> A_, final int LDA, final int TB, final int LTB, final Array<int> IPIV_, final int IPIV2, final Matrix<double> B_, final int LDB, final Array<double> WORK_, final int LWORK, final Box<int> INFO,) {
  final A = A_.dim();
  final IPIV = IPIV_.dim();
  final B = B_.dim();
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      String             UPLO;
      int                N, NRHS, LDA, LDB, LTB, LWORK, INFO;
      int                IPIV( * ), IPIV2( * );
      double             A( LDA, * ), B( LDB, * ), TB( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER, TQUERY, WQUERY;
      int                LWKMIN, LWKOPT;
r      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSYTRF_AA_2STAGE, DSYTRS_AA_2STAGE, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      WQUERY = ( LWORK == -1 );
      TQUERY = ( LTB == -1 );
      LWKMIN = max( 1, N );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LTB < max( 1, 4*N ) && !TQUERY ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -11;
      } else if ( LWORK < LWKMIN && !WQUERY ) {
         INFO = -13;
      }

      if ( INFO == 0 ) {
         dsytrf_aa_2stage(UPLO, N, A, LDA, TB, -1, IPIV, IPIV2, WORK, -1, INFO );
         LWKOPT = max( LWKMIN, INT( WORK( 1 ) ) );
         WORK[1] = LWKOPT;
      }

      if ( INFO != 0 ) {
         xerbla('DSYSV_AA_2STAGE', -INFO );
         return;
      } else if ( WQUERY || TQUERY ) {
         return;
      }

      // Compute the factorization A = U**T*T*U or A = L*T*L**T.

      dsytrf_aa_2stage(UPLO, N, A, LDA, TB, LTB, IPIV, IPIV2, WORK, LWORK, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         dsytrs_aa_2stage(UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, IPIV2, B, LDB, INFO );

      }

      WORK[1] = LWKOPT;

      }
