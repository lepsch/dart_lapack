import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgtsvx(final int FACT, final int TRANS, final int N, final int NRHS, final int DL, final int D, final int DU, final int DLF, final int DF, final int DUF, final int DU2, final Array<int> IPIV_, final Matrix<double> B_, final int LDB, final Matrix<double> X_, final int LDX, final int RCOND, final int FERR, final int BERR, final Array<double> _WORK_, final Array<int> IWORK_, final Box<int> INFO,) {
  final IPIV = IPIV_.dim();
  final B = B_.dim();
  final X = X_.dim();
  final _WORK = _WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             FACT, TRANS;
      int                INFO, LDB, LDX, N, NRHS;
      double             RCOND;
      int                IPIV( * ), IWORK( * );
      double             B( LDB, * ), BERR( * ), D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), FERR( * ), WORK( * ), X( LDX, * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      bool               NOFACT, NOTRAN;
      String             NORM;
      double             ANORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLANGT;
      // EXTERNAL lsame, DLAMCH, DLANGT
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGTCON, DGTRFS, DGTTRF, DGTTRS, DLACPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      NOTRAN = lsame( TRANS, 'N' );
      if ( !NOFACT && !lsame( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -14;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -16;
      }
      if ( INFO != 0 ) {
         xerbla('DGTSVX', -INFO );
         return;
      }

      if ( NOFACT ) {

         // Compute the LU factorization of A.

         dcopy(N, D, 1, DF, 1 );
         if ( N > 1 ) {
            dcopy(N-1, DL, 1, DLF, 1 );
            dcopy(N-1, DU, 1, DUF, 1 );
         }
         dgttrf(N, DLF, DF, DUF, DU2, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO;
            return;
         }
      }

      // Compute the norm of the matrix A.

      if ( NOTRAN ) {
         NORM = '1';
      } else {
         NORM = 'I';
      }
      ANORM = DLANGT( NORM, N, DL, D, DU );

      // Compute the reciprocal of the condition number of A.

      dgtcon(NORM, N, DLF, DF, DUF, DU2, IPIV, ANORM, RCOND, WORK, IWORK, INFO );

      // Compute the solution vectors X.

      dlacpy('Full', N, NRHS, B, LDB, X, LDX );
      dgttrs(TRANS, N, NRHS, DLF, DF, DUF, DU2, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      dgtrfs(TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO );

      // Set INFO = N+1 if the matrix is singular to working precision.

      if( RCOND < dlamch( 'Epsilon' ) ) INFO = N + 1;

      }
