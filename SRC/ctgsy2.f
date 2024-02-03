      SUBROUTINE CTGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N;
      REAL               RDSCAL, RDSUM, SCALE
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * ), D( LDD, * ), E( LDE, * ), F( LDF, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      int                LDZ;
      const              ZERO = 0.0, ONE = 1.0, LDZ = 2 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      int                I, IERR, J, K;
      REAL               SCALOC
      COMPLEX            ALPHA
      // ..
      // .. Local Arrays ..
      int                IPIV( LDZ ), JPIV( LDZ );
      COMPLEX            RHS( LDZ ), Z( LDZ, LDZ )
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CGESC2, CGETC2, CSCAL, CLATDF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, CONJG, MAX
      // ..
      // .. Executable Statements ..

      // Decode and test input parameters

      INFO = 0
      IERR = 0
      NOTRAN = LSAME( TRANS, 'N' )
      if ( .NOT.NOTRAN && .NOT.LSAME( TRANS, 'C' ) ) {
         INFO = -1
      } else if ( NOTRAN ) {
         if ( ( IJOB < 0 ) || ( IJOB > 2 ) ) {
            INFO = -2
         }
      }
      if ( INFO == 0 ) {
         if ( M <= 0 ) {
            INFO = -3
         } else if ( N <= 0 ) {
            INFO = -4
         } else if ( LDA < MAX( 1, M ) ) {
            INFO = -6
         } else if ( LDB < MAX( 1, N ) ) {
            INFO = -8
         } else if ( LDC < MAX( 1, M ) ) {
            INFO = -10
         } else if ( LDD < MAX( 1, M ) ) {
            INFO = -12
         } else if ( LDE < MAX( 1, N ) ) {
            INFO = -14
         } else if ( LDF < MAX( 1, M ) ) {
            INFO = -16
         }
      }
      if ( INFO != 0 ) {
         xerbla('CTGSY2', -INFO );
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

               cgetc2(LDZ, Z, LDZ, IPIV, JPIV, IERR );
               if (IERR > 0) INFO = IERR;
               if ( IJOB == 0 ) {
                  cgesc2(LDZ, Z, LDZ, RHS, IPIV, JPIV, SCALOC );
                  if ( SCALOC != ONE ) {
                     for (K = 1; K <= N; K++) { // 10
                        cscal(M, CMPLX( SCALOC, ZERO ), C( 1, K ), 1 );
                        cscal(M, CMPLX( SCALOC, ZERO ), F( 1, K ), 1 );
                     } // 10
                     SCALE = SCALE*SCALOC
                  }
               } else {
                  clatdf(IJOB, LDZ, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV );
               }

               // Unpack solution vector(s)

               C( I, J ) = RHS( 1 )
               F( I, J ) = RHS( 2 )

               // Substitute R(I, J) and L(I, J) into remaining equation.

               if ( I > 1 ) {
                  ALPHA = -RHS( 1 )
                  caxpy(I-1, ALPHA, A( 1, I ), 1, C( 1, J ), 1 );
                  caxpy(I-1, ALPHA, D( 1, I ), 1, F( 1, J ), 1 );
               }
               if ( J < N ) {
                  caxpy(N-J, RHS( 2 ), B( J, J+1 ), LDB, C( I, J+1 ), LDC );
                  caxpy(N-J, RHS( 2 ), E( J, J+1 ), LDE, F( I, J+1 ), LDF );
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

               Z( 1, 1 ) = CONJG( A( I, I ) )
               Z( 2, 1 ) = -CONJG( B( J, J ) )
               Z( 1, 2 ) = CONJG( D( I, I ) )
               Z( 2, 2 ) = -CONJG( E( J, J ) )


               // Set up right hand side(s)

               RHS( 1 ) = C( I, J )
               RHS( 2 ) = F( I, J )

               // Solve Z**H * x = RHS

               cgetc2(LDZ, Z, LDZ, IPIV, JPIV, IERR );
               if (IERR > 0) INFO = IERR;
               cgesc2(LDZ, Z, LDZ, RHS, IPIV, JPIV, SCALOC );
               if ( SCALOC != ONE ) {
                  for (K = 1; K <= N; K++) { // 40
                     cscal(M, CMPLX( SCALOC, ZERO ), C( 1, K ), 1 );
                     cscal(M, CMPLX( SCALOC, ZERO ), F( 1, K ), 1 );
                  } // 40
                  SCALE = SCALE*SCALOC
               }

               // Unpack solution vector(s)

               C( I, J ) = RHS( 1 )
               F( I, J ) = RHS( 2 )

               // Substitute R(I, J) and L(I, J) into remaining equation.

               for (K = 1; K <= J - 1; K++) { // 50
                  F( I, K ) = F( I, K ) + RHS( 1 )*CONJG( B( K, J ) ) + RHS( 2 )*CONJG( E( K, J ) )
               } // 50
               for (K = I + 1; K <= M; K++) { // 60
                  C( K, J ) = C( K, J ) - CONJG( A( I, K ) )*RHS( 1 ) - CONJG( D( I, K ) )*RHS( 2 )
               } // 60

            } // 70
         } // 80
      }
      RETURN

      // End of CTGSY2

      }
