import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dsysv(UPLO, N, NRHS, final Matrix<double> A, final int LDA, IPIV, final Matrix<double> B, final int LDB, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      int                IPIV( * );
      double             A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKOPT;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DSYTRF, DSYTRS, DSYTRS2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      if ( !lsame( UPLO, 'U' ) && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      } else if ( LWORK < 1 && !LQUERY ) {
         INFO = -10;
      }

      if ( INFO == 0 ) {
         if ( N == 0 ) {
            LWKOPT = 1;
         } else {
            dsytrf(UPLO, N, A, LDA, IPIV, WORK, -1, INFO );
            LWKOPT = INT( WORK( 1 ) );
         }
         WORK[1] = LWKOPT;
      }

      if ( INFO != 0 ) {
         xerbla('DSYSV ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Compute the factorization A = U*D*U**T or A = L*D*L**T.

      dsytrf(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         if ( LWORK < N ) {

         // Solve with TRS ( Use Level BLAS 2)

            dsytrs(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO );

         } else {

         // Solve with TRS2 ( Use Level BLAS 3)

            dsytrs2(UPLO,N,NRHS,A,LDA,IPIV,B,LDB,WORK,INFO );

         }

      }

      WORK[1] = LWKOPT;

      }
