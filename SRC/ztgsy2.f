      SUBROUTINE ZTGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N;
      double             RDSCAL, RDSUM, SCALE;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ), D( LDD, * ), E( LDE, * ), F( LDF, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      int                LDZ;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, LDZ = 2 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      int                I, IERR, J, K;
      double             SCALOC;
      COMPLEX*16         ALPHA
      // ..
      // .. Local Arrays ..
      int                IPIV( LDZ ), JPIV( LDZ );
      COMPLEX*16         RHS( LDZ ), Z( LDZ, LDZ )
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZGESC2, ZGETC2, ZLATDF, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, DCONJG, MAX
      // ..
      // .. Executable Statements ..

      // Decode and test input parameters

      INFO = 0
      IERR = 0
      NOTRAN = LSAME( TRANS, 'N' )
      if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) {
         INFO = -1
      } else if ( NOTRAN ) {
         if ( ( IJOB.LT.0 ) .OR. ( IJOB.GT.2 ) ) {
            INFO = -2
         }
      }
      if ( INFO == 0 ) {
         if ( M.LE.0 ) {
            INFO = -3
         } else if ( N.LE.0 ) {
            INFO = -4
         } else if ( LDA.LT.MAX( 1, M ) ) {
            INFO = -6
         } else if ( LDB.LT.MAX( 1, N ) ) {
            INFO = -8
         } else if ( LDC.LT.MAX( 1, M ) ) {
            INFO = -10
         } else if ( LDD.LT.MAX( 1, M ) ) {
            INFO = -12
         } else if ( LDE.LT.MAX( 1, N ) ) {
            INFO = -14
         } else if ( LDF.LT.MAX( 1, M ) ) {
            INFO = -16
         }
      }
      if ( INFO.NE.0 ) {
         xerbla('ZTGSY2', -INFO );
         RETURN
      }

      if ( NOTRAN ) {

         // Solve (I, J) - system
            // A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
            // D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
         // for I = M, M - 1, ..., 1; J = 1, 2, ..., N

         SCALE = ONE
         SCALOC = ONE
         for (J = 1; J <= N; J++) { // 30
            DO 20 I = M, 1, -1

               // Build 2 by 2 system

               Z( 1, 1 ) = A( I, I )
               Z( 2, 1 ) = D( I, I )
               Z( 1, 2 ) = -B( J, J )
               Z( 2, 2 ) = -E( J, J )

               // Set up right hand side(s)

               RHS( 1 ) = C( I, J )
               RHS( 2 ) = F( I, J )

               // Solve Z * x = RHS

               zgetc2(LDZ, Z, LDZ, IPIV, JPIV, IERR );
               if (IERR.GT.0) INFO = IERR;
               if ( IJOB == 0 ) {
                  zgesc2(LDZ, Z, LDZ, RHS, IPIV, JPIV, SCALOC );
                  if ( SCALOC.NE.ONE ) {
                     for (K = 1; K <= N; K++) { // 10
                        zscal(M, DCMPLX( SCALOC, ZERO ), C( 1, K ), 1 );
                        zscal(M, DCMPLX( SCALOC, ZERO ), F( 1, K ), 1 );
                     } // 10
                     SCALE = SCALE*SCALOC
                  }
               } else {
                  zlatdf(IJOB, LDZ, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV );
               }

               // Unpack solution vector(s)

               C( I, J ) = RHS( 1 )
               F( I, J ) = RHS( 2 )

               // Substitute R(I, J) and L(I, J) into remaining equation.

               if ( I.GT.1 ) {
                  ALPHA = -RHS( 1 )
                  zaxpy(I-1, ALPHA, A( 1, I ), 1, C( 1, J ), 1 );
                  zaxpy(I-1, ALPHA, D( 1, I ), 1, F( 1, J ), 1 );
               }
               if ( J.LT.N ) {
                  zaxpy(N-J, RHS( 2 ), B( J, J+1 ), LDB, C( I, J+1 ), LDC );
                  zaxpy(N-J, RHS( 2 ), E( J, J+1 ), LDE, F( I, J+1 ), LDF );
               }

            } // 20
         } // 30
      } else {

         // Solve transposed (I, J) - system:
            // A(I, I)**H * R(I, J) + D(I, I)**H * L(J, J) = C(I, J)
            // R(I, I) * B(J, J) + L(I, J) * E(J, J)   = -F(I, J)
         // for I = 1, 2, ..., M, J = N, N - 1, ..., 1

         SCALE = ONE
         SCALOC = ONE
         for (I = 1; I <= M; I++) { // 80
            DO 70 J = N, 1, -1

               // Build 2 by 2 system Z**H

               Z( 1, 1 ) = DCONJG( A( I, I ) )
               Z( 2, 1 ) = -DCONJG( B( J, J ) )
               Z( 1, 2 ) = DCONJG( D( I, I ) )
               Z( 2, 2 ) = -DCONJG( E( J, J ) )


               // Set up right hand side(s)

               RHS( 1 ) = C( I, J )
               RHS( 2 ) = F( I, J )

               // Solve Z**H * x = RHS

               zgetc2(LDZ, Z, LDZ, IPIV, JPIV, IERR );
               if (IERR.GT.0) INFO = IERR;
               zgesc2(LDZ, Z, LDZ, RHS, IPIV, JPIV, SCALOC );
               if ( SCALOC.NE.ONE ) {
                  for (K = 1; K <= N; K++) { // 40
                     zscal(M, DCMPLX( SCALOC, ZERO ), C( 1, K ), 1 );
                     zscal(M, DCMPLX( SCALOC, ZERO ), F( 1, K ), 1 );
                  } // 40
                  SCALE = SCALE*SCALOC
               }

               // Unpack solution vector(s)

               C( I, J ) = RHS( 1 )
               F( I, J ) = RHS( 2 )

               // Substitute R(I, J) and L(I, J) into remaining equation.

               for (K = 1; K <= J - 1; K++) { // 50
                  F( I, K ) = F( I, K ) + RHS( 1 )*DCONJG( B( K, J ) ) + RHS( 2 )*DCONJG( E( K, J ) )
               } // 50
               for (K = I + 1; K <= M; K++) { // 60
                  C( K, J ) = C( K, J ) - DCONJG( A( I, K ) )*RHS( 1 ) - DCONJG( D( I, K ) )*RHS( 2 )
               } // 60

            } // 70
         } // 80
      }
      RETURN

      // End of ZTGSY2

      }
