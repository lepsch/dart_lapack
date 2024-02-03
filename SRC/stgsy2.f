      SUBROUTINE STGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL, IWORK, PQ, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N, PQ;
      REAL               RDSCAL, RDSUM, SCALE
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * ), D( LDD, * ), E( LDE, * ), F( LDF, * )
      // ..

*  =====================================================================
*  Replaced various illegal calls to SCOPY by calls to SLASET.
*  Sven Hammarling, 27/5/02.

      // .. Parameters ..
      int                LDZ;
      const              LDZ = 8 ;
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      int                I, IE, IERR, II, IS, ISP1, J, JE, JJ, JS, JSP1, K, MB, NB, P, Q, ZDIM;
      REAL               ALPHA, SCALOC
      // ..
      // .. Local Arrays ..
      int                IPIV( LDZ ), JPIV( LDZ );
      REAL               RHS( LDZ ), Z( LDZ, LDZ )
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SGEMM, SGEMV, SGER, SGESC2, SGETC2, SSCAL, SLASET, SLATDF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Decode and test input parameters

      INFO = 0
      IERR = 0
      NOTRAN = LSAME( TRANS, 'N' )
      if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) {
         INFO = -1
      } else if ( NOTRAN ) {
         if ( ( IJOB.LT.0 ) .OR. ( IJOB.GT.2 ) ) {
            INFO = -2
         }
      }
      if ( INFO.EQ.0 ) {
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
         xerbla('STGSY2', -INFO );
         RETURN
      }

      // Determine block structure of A

      PQ = 0
      P = 0
      I = 1
      } // 10
      IF( I.GT.M ) GO TO 20
      P = P + 1
      IWORK( P ) = I
      IF( I.EQ.M ) GO TO 20
      if ( A( I+1, I ).NE.ZERO ) {
         I = I + 2
      } else {
         I = I + 1
      }
      GO TO 10
      } // 20
      IWORK( P+1 ) = M + 1

      // Determine block structure of B

      Q = P + 1
      J = 1
      } // 30
      IF( J.GT.N ) GO TO 40
      Q = Q + 1
      IWORK( Q ) = J
      IF( J.EQ.N ) GO TO 40
      if ( B( J+1, J ).NE.ZERO ) {
         J = J + 2
      } else {
         J = J + 1
      }
      GO TO 30
      } // 40
      IWORK( Q+1 ) = N + 1
      PQ = P*( Q-P-1 )

      if ( NOTRAN ) {

         // Solve (I, J) - subsystem
            // A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
            // D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
         // for I = P, P - 1, ..., 1; J = 1, 2, ..., Q

         SCALE = ONE
         SCALOC = ONE
         for (J = P + 2; J <= Q; J++) { // 120
            JS = IWORK( J )
            JSP1 = JS + 1
            JE = IWORK( J+1 ) - 1
            NB = JE - JS + 1
            DO 110 I = P, 1, -1

               IS = IWORK( I )
               ISP1 = IS + 1
               IE = IWORK( I+1 ) - 1
               MB = IE - IS + 1
               ZDIM = MB*NB*2

               if ( ( MB.EQ.1 ) .AND. ( NB.EQ.1 ) ) {

                  // Build a 2-by-2 system Z * x = RHS

                  Z( 1, 1 ) = A( IS, IS )
                  Z( 2, 1 ) = D( IS, IS )
                  Z( 1, 2 ) = -B( JS, JS )
                  Z( 2, 2 ) = -E( JS, JS )

                  // Set up right hand side(s)

                  RHS( 1 ) = C( IS, JS )
                  RHS( 2 ) = F( IS, JS )

                  // Solve Z * x = RHS

                  sgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR );
                  IF( IERR.GT.0 ) INFO = IERR

                  if ( IJOB.EQ.0 ) {
                     sgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC );
                     if ( SCALOC.NE.ONE ) {
                        for (K = 1; K <= N; K++) { // 50
                           sscal(M, SCALOC, C( 1, K ), 1 );
                           sscal(M, SCALOC, F( 1, K ), 1 );
                        } // 50
                        SCALE = SCALE*SCALOC
                     }
                  } else {
                     slatdf(IJOB, ZDIM, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV );
                  }

                  // Unpack solution vector(s)

                  C( IS, JS ) = RHS( 1 )
                  F( IS, JS ) = RHS( 2 )

                  // Substitute R(I, J) and L(I, J) into remaining
                  // equation.

                  if ( I.GT.1 ) {
                     ALPHA = -RHS( 1 )
                     saxpy(IS-1, ALPHA, A( 1, IS ), 1, C( 1, JS ), 1 )                      CALL SAXPY( IS-1, ALPHA, D( 1, IS ), 1, F( 1, JS ), 1 );
                  }
                  if ( J.LT.Q ) {
                     saxpy(N-JE, RHS( 2 ), B( JS, JE+1 ), LDB, C( IS, JE+1 ), LDC )                      CALL SAXPY( N-JE, RHS( 2 ), E( JS, JE+1 ), LDE, F( IS, JE+1 ), LDF );
                  }

               } else if ( ( MB.EQ.1 ) .AND. ( NB.EQ.2 ) ) {

                  // Build a 4-by-4 system Z * x = RHS

                  Z( 1, 1 ) = A( IS, IS )
                  Z( 2, 1 ) = ZERO
                  Z( 3, 1 ) = D( IS, IS )
                  Z( 4, 1 ) = ZERO

                  Z( 1, 2 ) = ZERO
                  Z( 2, 2 ) = A( IS, IS )
                  Z( 3, 2 ) = ZERO
                  Z( 4, 2 ) = D( IS, IS )

                  Z( 1, 3 ) = -B( JS, JS )
                  Z( 2, 3 ) = -B( JS, JSP1 )
                  Z( 3, 3 ) = -E( JS, JS )
                  Z( 4, 3 ) = -E( JS, JSP1 )

                  Z( 1, 4 ) = -B( JSP1, JS )
                  Z( 2, 4 ) = -B( JSP1, JSP1 )
                  Z( 3, 4 ) = ZERO
                  Z( 4, 4 ) = -E( JSP1, JSP1 )

                  // Set up right hand side(s)

                  RHS( 1 ) = C( IS, JS )
                  RHS( 2 ) = C( IS, JSP1 )
                  RHS( 3 ) = F( IS, JS )
                  RHS( 4 ) = F( IS, JSP1 )

                  // Solve Z * x = RHS

                  sgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR );
                  IF( IERR.GT.0 ) INFO = IERR

                  if ( IJOB.EQ.0 ) {
                     sgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC );
                     if ( SCALOC.NE.ONE ) {
                        for (K = 1; K <= N; K++) { // 60
                           sscal(M, SCALOC, C( 1, K ), 1 );
                           sscal(M, SCALOC, F( 1, K ), 1 );
                        } // 60
                        SCALE = SCALE*SCALOC
                     }
                  } else {
                     slatdf(IJOB, ZDIM, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV );
                  }

                  // Unpack solution vector(s)

                  C( IS, JS ) = RHS( 1 )
                  C( IS, JSP1 ) = RHS( 2 )
                  F( IS, JS ) = RHS( 3 )
                  F( IS, JSP1 ) = RHS( 4 )

                  // Substitute R(I, J) and L(I, J) into remaining
                  // equation.

                  if ( I.GT.1 ) {
                     sger(IS-1, NB, -ONE, A( 1, IS ), 1, RHS( 1 ), 1, C( 1, JS ), LDC )                      CALL SGER( IS-1, NB, -ONE, D( 1, IS ), 1, RHS( 1 ), 1, F( 1, JS ), LDF );
                  }
                  if ( J.LT.Q ) {
                     saxpy(N-JE, RHS( 3 ), B( JS, JE+1 ), LDB, C( IS, JE+1 ), LDC )                      CALL SAXPY( N-JE, RHS( 3 ), E( JS, JE+1 ), LDE, F( IS, JE+1 ), LDF )                      CALL SAXPY( N-JE, RHS( 4 ), B( JSP1, JE+1 ), LDB, C( IS, JE+1 ), LDC )                      CALL SAXPY( N-JE, RHS( 4 ), E( JSP1, JE+1 ), LDE, F( IS, JE+1 ), LDF );
                  }

               } else if ( ( MB.EQ.2 ) .AND. ( NB.EQ.1 ) ) {

                  // Build a 4-by-4 system Z * x = RHS

                  Z( 1, 1 ) = A( IS, IS )
                  Z( 2, 1 ) = A( ISP1, IS )
                  Z( 3, 1 ) = D( IS, IS )
                  Z( 4, 1 ) = ZERO

                  Z( 1, 2 ) = A( IS, ISP1 )
                  Z( 2, 2 ) = A( ISP1, ISP1 )
                  Z( 3, 2 ) = D( IS, ISP1 )
                  Z( 4, 2 ) = D( ISP1, ISP1 )

                  Z( 1, 3 ) = -B( JS, JS )
                  Z( 2, 3 ) = ZERO
                  Z( 3, 3 ) = -E( JS, JS )
                  Z( 4, 3 ) = ZERO

                  Z( 1, 4 ) = ZERO
                  Z( 2, 4 ) = -B( JS, JS )
                  Z( 3, 4 ) = ZERO
                  Z( 4, 4 ) = -E( JS, JS )

                  // Set up right hand side(s)

                  RHS( 1 ) = C( IS, JS )
                  RHS( 2 ) = C( ISP1, JS )
                  RHS( 3 ) = F( IS, JS )
                  RHS( 4 ) = F( ISP1, JS )

                  // Solve Z * x = RHS

                  sgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR );
                  IF( IERR.GT.0 ) INFO = IERR
                  if ( IJOB.EQ.0 ) {
                     sgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC );
                     if ( SCALOC.NE.ONE ) {
                        for (K = 1; K <= N; K++) { // 70
                           sscal(M, SCALOC, C( 1, K ), 1 );
                           sscal(M, SCALOC, F( 1, K ), 1 );
                        } // 70
                        SCALE = SCALE*SCALOC
                     }
                  } else {
                     slatdf(IJOB, ZDIM, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV );
                  }

                  // Unpack solution vector(s)

                  C( IS, JS ) = RHS( 1 )
                  C( ISP1, JS ) = RHS( 2 )
                  F( IS, JS ) = RHS( 3 )
                  F( ISP1, JS ) = RHS( 4 )

                  // Substitute R(I, J) and L(I, J) into remaining
                  // equation.

                  if ( I.GT.1 ) {
                     sgemv('N', IS-1, MB, -ONE, A( 1, IS ), LDA, RHS( 1 ), 1, ONE, C( 1, JS ), 1 )                      CALL SGEMV( 'N', IS-1, MB, -ONE, D( 1, IS ), LDD, RHS( 1 ), 1, ONE, F( 1, JS ), 1 );
                  }
                  if ( J.LT.Q ) {
                     sger(MB, N-JE, ONE, RHS( 3 ), 1, B( JS, JE+1 ), LDB, C( IS, JE+1 ), LDC )                      CALL SGER( MB, N-JE, ONE, RHS( 3 ), 1, E( JS, JE+1 ), LDE, F( IS, JE+1 ), LDF );
                  }

               } else if ( ( MB.EQ.2 ) .AND. ( NB.EQ.2 ) ) {

                  // Build an 8-by-8 system Z * x = RHS

                  slaset('F', LDZ, LDZ, ZERO, ZERO, Z, LDZ );

                  Z( 1, 1 ) = A( IS, IS )
                  Z( 2, 1 ) = A( ISP1, IS )
                  Z( 5, 1 ) = D( IS, IS )

                  Z( 1, 2 ) = A( IS, ISP1 )
                  Z( 2, 2 ) = A( ISP1, ISP1 )
                  Z( 5, 2 ) = D( IS, ISP1 )
                  Z( 6, 2 ) = D( ISP1, ISP1 )

                  Z( 3, 3 ) = A( IS, IS )
                  Z( 4, 3 ) = A( ISP1, IS )
                  Z( 7, 3 ) = D( IS, IS )

                  Z( 3, 4 ) = A( IS, ISP1 )
                  Z( 4, 4 ) = A( ISP1, ISP1 )
                  Z( 7, 4 ) = D( IS, ISP1 )
                  Z( 8, 4 ) = D( ISP1, ISP1 )

                  Z( 1, 5 ) = -B( JS, JS )
                  Z( 3, 5 ) = -B( JS, JSP1 )
                  Z( 5, 5 ) = -E( JS, JS )
                  Z( 7, 5 ) = -E( JS, JSP1 )

                  Z( 2, 6 ) = -B( JS, JS )
                  Z( 4, 6 ) = -B( JS, JSP1 )
                  Z( 6, 6 ) = -E( JS, JS )
                  Z( 8, 6 ) = -E( JS, JSP1 )

                  Z( 1, 7 ) = -B( JSP1, JS )
                  Z( 3, 7 ) = -B( JSP1, JSP1 )
                  Z( 7, 7 ) = -E( JSP1, JSP1 )

                  Z( 2, 8 ) = -B( JSP1, JS )
                  Z( 4, 8 ) = -B( JSP1, JSP1 )
                  Z( 8, 8 ) = -E( JSP1, JSP1 )

                  // Set up right hand side(s)

                  K = 1
                  II = MB*NB + 1
                  for (JJ = 0; JJ <= NB - 1; JJ++) { // 80
                     scopy(MB, C( IS, JS+JJ ), 1, RHS( K ), 1 );
                     scopy(MB, F( IS, JS+JJ ), 1, RHS( II ), 1 );
                     K = K + MB
                     II = II + MB
                  } // 80

                  // Solve Z * x = RHS

                  sgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR );
                  IF( IERR.GT.0 ) INFO = IERR
                  if ( IJOB.EQ.0 ) {
                     sgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC );
                     if ( SCALOC.NE.ONE ) {
                        for (K = 1; K <= N; K++) { // 90
                           sscal(M, SCALOC, C( 1, K ), 1 );
                           sscal(M, SCALOC, F( 1, K ), 1 );
                        } // 90
                        SCALE = SCALE*SCALOC
                     }
                  } else {
                     slatdf(IJOB, ZDIM, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV );
                  }

                  // Unpack solution vector(s)

                  K = 1
                  II = MB*NB + 1
                  for (JJ = 0; JJ <= NB - 1; JJ++) { // 100
                     scopy(MB, RHS( K ), 1, C( IS, JS+JJ ), 1 );
                     scopy(MB, RHS( II ), 1, F( IS, JS+JJ ), 1 );
                     K = K + MB
                     II = II + MB
                  } // 100

                  // Substitute R(I, J) and L(I, J) into remaining
                  // equation.

                  if ( I.GT.1 ) {
                     sgemm('N', 'N', IS-1, NB, MB, -ONE, A( 1, IS ), LDA, RHS( 1 ), MB, ONE, C( 1, JS ), LDC )                      CALL SGEMM( 'N', 'N', IS-1, NB, MB, -ONE, D( 1, IS ), LDD, RHS( 1 ), MB, ONE, F( 1, JS ), LDF );
                  }
                  if ( J.LT.Q ) {
                     K = MB*NB + 1
                     sgemm('N', 'N', MB, N-JE, NB, ONE, RHS( K ), MB, B( JS, JE+1 ), LDB, ONE, C( IS, JE+1 ), LDC )                      CALL SGEMM( 'N', 'N', MB, N-JE, NB, ONE, RHS( K ), MB, E( JS, JE+1 ), LDE, ONE, F( IS, JE+1 ), LDF );
                  }

               }

            } // 110
         } // 120
      } else {

         // Solve (I, J) - subsystem
              // A(I, I)**T * R(I, J) + D(I, I)**T * L(J, J)  =  C(I, J)
              // R(I, I)  * B(J, J) + L(I, J)  * E(J, J)  = -F(I, J)
         // for I = 1, 2, ..., P, J = Q, Q - 1, ..., 1

         SCALE = ONE
         SCALOC = ONE
         for (I = 1; I <= P; I++) { // 200

            IS = IWORK( I )
            ISP1 = IS + 1
            IE = IWORK( I+1 ) - 1
            MB = IE - IS + 1
            DO 190 J = Q, P + 2, -1

               JS = IWORK( J )
               JSP1 = JS + 1
               JE = IWORK( J+1 ) - 1
               NB = JE - JS + 1
               ZDIM = MB*NB*2
               if ( ( MB.EQ.1 ) .AND. ( NB.EQ.1 ) ) {

                  // Build a 2-by-2 system Z**T * x = RHS

                  Z( 1, 1 ) = A( IS, IS )
                  Z( 2, 1 ) = -B( JS, JS )
                  Z( 1, 2 ) = D( IS, IS )
                  Z( 2, 2 ) = -E( JS, JS )

                  // Set up right hand side(s)

                  RHS( 1 ) = C( IS, JS )
                  RHS( 2 ) = F( IS, JS )

                  // Solve Z**T * x = RHS

                  sgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR );
                  IF( IERR.GT.0 ) INFO = IERR

                  sgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC );
                  if ( SCALOC.NE.ONE ) {
                     for (K = 1; K <= N; K++) { // 130
                        sscal(M, SCALOC, C( 1, K ), 1 );
                        sscal(M, SCALOC, F( 1, K ), 1 );
                     } // 130
                     SCALE = SCALE*SCALOC
                  }

                  // Unpack solution vector(s)

                  C( IS, JS ) = RHS( 1 )
                  F( IS, JS ) = RHS( 2 )

                  // Substitute R(I, J) and L(I, J) into remaining
                  // equation.

                  if ( J.GT.P+2 ) {
                     ALPHA = RHS( 1 )
                     saxpy(JS-1, ALPHA, B( 1, JS ), 1, F( IS, 1 ), LDF );
                     ALPHA = RHS( 2 )
                     saxpy(JS-1, ALPHA, E( 1, JS ), 1, F( IS, 1 ), LDF );
                  }
                  if ( I.LT.P ) {
                     ALPHA = -RHS( 1 )
                     saxpy(M-IE, ALPHA, A( IS, IE+1 ), LDA, C( IE+1, JS ), 1 );
                     ALPHA = -RHS( 2 )
                     saxpy(M-IE, ALPHA, D( IS, IE+1 ), LDD, C( IE+1, JS ), 1 );
                  }

               } else if ( ( MB.EQ.1 ) .AND. ( NB.EQ.2 ) ) {

                  // Build a 4-by-4 system Z**T * x = RHS

                  Z( 1, 1 ) = A( IS, IS )
                  Z( 2, 1 ) = ZERO
                  Z( 3, 1 ) = -B( JS, JS )
                  Z( 4, 1 ) = -B( JSP1, JS )

                  Z( 1, 2 ) = ZERO
                  Z( 2, 2 ) = A( IS, IS )
                  Z( 3, 2 ) = -B( JS, JSP1 )
                  Z( 4, 2 ) = -B( JSP1, JSP1 )

                  Z( 1, 3 ) = D( IS, IS )
                  Z( 2, 3 ) = ZERO
                  Z( 3, 3 ) = -E( JS, JS )
                  Z( 4, 3 ) = ZERO

                  Z( 1, 4 ) = ZERO
                  Z( 2, 4 ) = D( IS, IS )
                  Z( 3, 4 ) = -E( JS, JSP1 )
                  Z( 4, 4 ) = -E( JSP1, JSP1 )

                  // Set up right hand side(s)

                  RHS( 1 ) = C( IS, JS )
                  RHS( 2 ) = C( IS, JSP1 )
                  RHS( 3 ) = F( IS, JS )
                  RHS( 4 ) = F( IS, JSP1 )

                  // Solve Z**T * x = RHS

                  sgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR );
                  IF( IERR.GT.0 ) INFO = IERR
                  sgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC );
                  if ( SCALOC.NE.ONE ) {
                     for (K = 1; K <= N; K++) { // 140
                        sscal(M, SCALOC, C( 1, K ), 1 );
                        sscal(M, SCALOC, F( 1, K ), 1 );
                     } // 140
                     SCALE = SCALE*SCALOC
                  }

                  // Unpack solution vector(s)

                  C( IS, JS ) = RHS( 1 )
                  C( IS, JSP1 ) = RHS( 2 )
                  F( IS, JS ) = RHS( 3 )
                  F( IS, JSP1 ) = RHS( 4 )

                  // Substitute R(I, J) and L(I, J) into remaining
                  // equation.

                  if ( J.GT.P+2 ) {
                     saxpy(JS-1, RHS( 1 ), B( 1, JS ), 1, F( IS, 1 ), LDF )                      CALL SAXPY( JS-1, RHS( 2 ), B( 1, JSP1 ), 1, F( IS, 1 ), LDF )                      CALL SAXPY( JS-1, RHS( 3 ), E( 1, JS ), 1, F( IS, 1 ), LDF )                      CALL SAXPY( JS-1, RHS( 4 ), E( 1, JSP1 ), 1, F( IS, 1 ), LDF );
                  }
                  if ( I.LT.P ) {
                     sger(M-IE, NB, -ONE, A( IS, IE+1 ), LDA, RHS( 1 ), 1, C( IE+1, JS ), LDC )                      CALL SGER( M-IE, NB, -ONE, D( IS, IE+1 ), LDD, RHS( 3 ), 1, C( IE+1, JS ), LDC );
                  }

               } else if ( ( MB.EQ.2 ) .AND. ( NB.EQ.1 ) ) {

                  // Build a 4-by-4 system Z**T * x = RHS

                  Z( 1, 1 ) = A( IS, IS )
                  Z( 2, 1 ) = A( IS, ISP1 )
                  Z( 3, 1 ) = -B( JS, JS )
                  Z( 4, 1 ) = ZERO

                  Z( 1, 2 ) = A( ISP1, IS )
                  Z( 2, 2 ) = A( ISP1, ISP1 )
                  Z( 3, 2 ) = ZERO
                  Z( 4, 2 ) = -B( JS, JS )

                  Z( 1, 3 ) = D( IS, IS )
                  Z( 2, 3 ) = D( IS, ISP1 )
                  Z( 3, 3 ) = -E( JS, JS )
                  Z( 4, 3 ) = ZERO

                  Z( 1, 4 ) = ZERO
                  Z( 2, 4 ) = D( ISP1, ISP1 )
                  Z( 3, 4 ) = ZERO
                  Z( 4, 4 ) = -E( JS, JS )

                  // Set up right hand side(s)

                  RHS( 1 ) = C( IS, JS )
                  RHS( 2 ) = C( ISP1, JS )
                  RHS( 3 ) = F( IS, JS )
                  RHS( 4 ) = F( ISP1, JS )

                  // Solve Z**T * x = RHS

                  sgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR );
                  IF( IERR.GT.0 ) INFO = IERR

                  sgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC );
                  if ( SCALOC.NE.ONE ) {
                     for (K = 1; K <= N; K++) { // 150
                        sscal(M, SCALOC, C( 1, K ), 1 );
                        sscal(M, SCALOC, F( 1, K ), 1 );
                     } // 150
                     SCALE = SCALE*SCALOC
                  }

                  // Unpack solution vector(s)

                  C( IS, JS ) = RHS( 1 )
                  C( ISP1, JS ) = RHS( 2 )
                  F( IS, JS ) = RHS( 3 )
                  F( ISP1, JS ) = RHS( 4 )

                  // Substitute R(I, J) and L(I, J) into remaining
                  // equation.

                  if ( J.GT.P+2 ) {
                     sger(MB, JS-1, ONE, RHS( 1 ), 1, B( 1, JS ), 1, F( IS, 1 ), LDF )                      CALL SGER( MB, JS-1, ONE, RHS( 3 ), 1, E( 1, JS ), 1, F( IS, 1 ), LDF );
                  }
                  if ( I.LT.P ) {
                     sgemv('T', MB, M-IE, -ONE, A( IS, IE+1 ), LDA, RHS( 1 ), 1, ONE, C( IE+1, JS ), 1 )                      CALL SGEMV( 'T', MB, M-IE, -ONE, D( IS, IE+1 ), LDD, RHS( 3 ), 1, ONE, C( IE+1, JS ), 1 );
                  }

               } else if ( ( MB.EQ.2 ) .AND. ( NB.EQ.2 ) ) {

                  // Build an 8-by-8 system Z**T * x = RHS

                  slaset('F', LDZ, LDZ, ZERO, ZERO, Z, LDZ );

                  Z( 1, 1 ) = A( IS, IS )
                  Z( 2, 1 ) = A( IS, ISP1 )
                  Z( 5, 1 ) = -B( JS, JS )
                  Z( 7, 1 ) = -B( JSP1, JS )

                  Z( 1, 2 ) = A( ISP1, IS )
                  Z( 2, 2 ) = A( ISP1, ISP1 )
                  Z( 6, 2 ) = -B( JS, JS )
                  Z( 8, 2 ) = -B( JSP1, JS )

                  Z( 3, 3 ) = A( IS, IS )
                  Z( 4, 3 ) = A( IS, ISP1 )
                  Z( 5, 3 ) = -B( JS, JSP1 )
                  Z( 7, 3 ) = -B( JSP1, JSP1 )

                  Z( 3, 4 ) = A( ISP1, IS )
                  Z( 4, 4 ) = A( ISP1, ISP1 )
                  Z( 6, 4 ) = -B( JS, JSP1 )
                  Z( 8, 4 ) = -B( JSP1, JSP1 )

                  Z( 1, 5 ) = D( IS, IS )
                  Z( 2, 5 ) = D( IS, ISP1 )
                  Z( 5, 5 ) = -E( JS, JS )

                  Z( 2, 6 ) = D( ISP1, ISP1 )
                  Z( 6, 6 ) = -E( JS, JS )

                  Z( 3, 7 ) = D( IS, IS )
                  Z( 4, 7 ) = D( IS, ISP1 )
                  Z( 5, 7 ) = -E( JS, JSP1 )
                  Z( 7, 7 ) = -E( JSP1, JSP1 )

                  Z( 4, 8 ) = D( ISP1, ISP1 )
                  Z( 6, 8 ) = -E( JS, JSP1 )
                  Z( 8, 8 ) = -E( JSP1, JSP1 )

                  // Set up right hand side(s)

                  K = 1
                  II = MB*NB + 1
                  for (JJ = 0; JJ <= NB - 1; JJ++) { // 160
                     scopy(MB, C( IS, JS+JJ ), 1, RHS( K ), 1 );
                     scopy(MB, F( IS, JS+JJ ), 1, RHS( II ), 1 );
                     K = K + MB
                     II = II + MB
                  } // 160


                  // Solve Z**T * x = RHS

                  sgetc2(ZDIM, Z, LDZ, IPIV, JPIV, IERR );
                  IF( IERR.GT.0 ) INFO = IERR

                  sgesc2(ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC );
                  if ( SCALOC.NE.ONE ) {
                     for (K = 1; K <= N; K++) { // 170
                        sscal(M, SCALOC, C( 1, K ), 1 );
                        sscal(M, SCALOC, F( 1, K ), 1 );
                     } // 170
                     SCALE = SCALE*SCALOC
                  }

                  // Unpack solution vector(s)

                  K = 1
                  II = MB*NB + 1
                  for (JJ = 0; JJ <= NB - 1; JJ++) { // 180
                     scopy(MB, RHS( K ), 1, C( IS, JS+JJ ), 1 );
                     scopy(MB, RHS( II ), 1, F( IS, JS+JJ ), 1 );
                     K = K + MB
                     II = II + MB
                  } // 180

                  // Substitute R(I, J) and L(I, J) into remaining
                  // equation.

                  if ( J.GT.P+2 ) {
                     sgemm('N', 'T', MB, JS-1, NB, ONE, C( IS, JS ), LDC, B( 1, JS ), LDB, ONE, F( IS, 1 ), LDF )                      CALL SGEMM( 'N', 'T', MB, JS-1, NB, ONE, F( IS, JS ), LDF, E( 1, JS ), LDE, ONE, F( IS, 1 ), LDF );
                  }
                  if ( I.LT.P ) {
                     sgemm('T', 'N', M-IE, NB, MB, -ONE, A( IS, IE+1 ), LDA, C( IS, JS ), LDC, ONE, C( IE+1, JS ), LDC )                      CALL SGEMM( 'T', 'N', M-IE, NB, MB, -ONE, D( IS, IE+1 ), LDD, F( IS, JS ), LDF, ONE, C( IE+1, JS ), LDC );
                  }

               }

            } // 190
         } // 200

      }
      RETURN

      // End of STGSY2

      }
