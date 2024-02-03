      SUBROUTINE CLAROT( LROWS, LLEFT, LRIGHT, NL, C, S, A, LDA, XLEFT, XRIGHT )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LLEFT, LRIGHT, LROWS;
      int                LDA, NL;
      COMPLEX            C, S, XLEFT, XRIGHT
      // ..
      // .. Array Arguments ..
      COMPLEX            A( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                IINC, INEXT, IX, IY, IYT, J, NT;
      COMPLEX            TEMPX
      // ..
      // .. Local Arrays ..
      COMPLEX            XT( 2 ), YT( 2 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG
      // ..
      // .. Executable Statements ..

      // Set up indices, arrays for ends

      IF( LROWS ) THEN
         IINC = LDA
         INEXT = 1
      ELSE
         IINC = 1
         INEXT = LDA
      END IF

      IF( LLEFT ) THEN
         NT = 1
         IX = 1 + IINC
         IY = 2 + LDA
         XT( 1 ) = A( 1 )
         YT( 1 ) = XLEFT
      ELSE
         NT = 0
         IX = 1
         IY = 1 + INEXT
      END IF

      IF( LRIGHT ) THEN
         IYT = 1 + INEXT + ( NL-1 )*IINC
         NT = NT + 1
         XT( NT ) = XRIGHT
         YT( NT ) = A( IYT )
      END IF

      // Check for errors

      IF( NL.LT.NT ) THEN
         CALL XERBLA( 'CLAROT', 4 )
         RETURN
      END IF
      IF( LDA.LE.0 .OR. ( .NOT.LROWS .AND. LDA.LT.NL-NT ) ) THEN
         CALL XERBLA( 'CLAROT', 8 )
         RETURN
      END IF

      // Rotate

      // CROT( NL-NT, A(IX),IINC, A(IY),IINC, C, S ) with complex C, S

      DO 10 J = 0, NL - NT - 1
         TEMPX = C*A( IX+J*IINC ) + S*A( IY+J*IINC )
         A( IY+J*IINC ) = -CONJG( S )*A( IX+J*IINC ) + CONJG( C )*A( IY+J*IINC )
         A( IX+J*IINC ) = TEMPX
   10 CONTINUE

      // CROT( NT, XT,1, YT,1, C, S ) with complex C, S

      DO 20 J = 1, NT
         TEMPX = C*XT( J ) + S*YT( J )
         YT( J ) = -CONJG( S )*XT( J ) + CONJG( C )*YT( J )
         XT( J ) = TEMPX
   20 CONTINUE

      // Stuff values back into XLEFT, XRIGHT, etc.

      IF( LLEFT ) THEN
         A( 1 ) = XT( 1 )
         XLEFT = YT( 1 )
      END IF

      IF( LRIGHT ) THEN
         XRIGHT = XT( NT )
         A( IYT ) = YT( NT )
      END IF

      RETURN

      // End of CLAROT

      END
