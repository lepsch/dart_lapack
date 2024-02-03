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
      const              ZERO = 0.0E+0, ONE = 1.0E+0, LDZ = 2 ;
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
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( NOTRAN ) THEN
         IF( ( IJOB.LT.0 ) .OR. ( IJOB.GT.2 ) ) THEN
            INFO = -2
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( M.LE.0 ) THEN
            INFO = -3
         ELSE IF( N.LE.0 ) THEN
            INFO = -4
         ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
            INFO = -6
         ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
            INFO = -8
         ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
            INFO = -10
         ELSE IF( LDD.LT.MAX( 1, M ) ) THEN
            INFO = -12
         ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
            INFO = -14
         ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
            INFO = -16
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTGSY2', -INFO )
         RETURN
      END IF

      IF( NOTRAN ) THEN

         // Solve (I, J) - system
            // A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
            // D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
         // for I = M, M - 1, ..., 1; J = 1, 2, ..., N

         SCALE = ONE
         SCALOC = ONE
         DO 30 J = 1, N
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

               CALL CGETC2( LDZ, Z, LDZ, IPIV, JPIV, IERR )
               IF( IERR.GT.0 ) INFO = IERR
               IF( IJOB.EQ.0 ) THEN
                  CALL CGESC2( LDZ, Z, LDZ, RHS, IPIV, JPIV, SCALOC )
                  IF( SCALOC.NE.ONE ) THEN
                     DO 10 K = 1, N
                        CALL CSCAL( M, CMPLX( SCALOC, ZERO ), C( 1, K ), 1 )                         CALL CSCAL( M, CMPLX( SCALOC, ZERO ), F( 1, K ), 1 )
   10                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
               ELSE
                  CALL CLATDF( IJOB, LDZ, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV )
               END IF

               // Unpack solution vector(s)

               C( I, J ) = RHS( 1 )
               F( I, J ) = RHS( 2 )

               // Substitute R(I, J) and L(I, J) into remaining equation.

               IF( I.GT.1 ) THEN
                  ALPHA = -RHS( 1 )
                  CALL CAXPY( I-1, ALPHA, A( 1, I ), 1, C( 1, J ), 1 )
                  CALL CAXPY( I-1, ALPHA, D( 1, I ), 1, F( 1, J ), 1 )
               END IF
               IF( J.LT.N ) THEN
                  CALL CAXPY( N-J, RHS( 2 ), B( J, J+1 ), LDB, C( I, J+1 ), LDC )                   CALL CAXPY( N-J, RHS( 2 ), E( J, J+1 ), LDE, F( I, J+1 ), LDF )
               END IF

   20       CONTINUE
   30    CONTINUE
      ELSE

         // Solve transposed (I, J) - system:
            // A(I, I)**H * R(I, J) + D(I, I)**H * L(J, J) = C(I, J)
            // R(I, I) * B(J, J) + L(I, J) * E(J, J)   = -F(I, J)
         // for I = 1, 2, ..., M, J = N, N - 1, ..., 1

         SCALE = ONE
         SCALOC = ONE
         DO 80 I = 1, M
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

               CALL CGETC2( LDZ, Z, LDZ, IPIV, JPIV, IERR )
               IF( IERR.GT.0 ) INFO = IERR
               CALL CGESC2( LDZ, Z, LDZ, RHS, IPIV, JPIV, SCALOC )
               IF( SCALOC.NE.ONE ) THEN
                  DO 40 K = 1, N
                     CALL CSCAL( M, CMPLX( SCALOC, ZERO ), C( 1, K ), 1 )                      CALL CSCAL( M, CMPLX( SCALOC, ZERO ), F( 1, K ), 1 )
   40             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF

               // Unpack solution vector(s)

               C( I, J ) = RHS( 1 )
               F( I, J ) = RHS( 2 )

               // Substitute R(I, J) and L(I, J) into remaining equation.

               DO 50 K = 1, J - 1
                  F( I, K ) = F( I, K ) + RHS( 1 )*CONJG( B( K, J ) ) + RHS( 2 )*CONJG( E( K, J ) )
   50          CONTINUE
               DO 60 K = I + 1, M
                  C( K, J ) = C( K, J ) - CONJG( A( I, K ) )*RHS( 1 ) - CONJG( D( I, K ) )*RHS( 2 )
   60          CONTINUE

   70       CONTINUE
   80    CONTINUE
      END IF
      RETURN

      // End of CTGSY2

      }
