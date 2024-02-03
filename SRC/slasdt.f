      SUBROUTINE SLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                LVL, MSUB, N, ND;
      // ..
      // .. Array Arguments ..
      int                INODE( * ), NDIML( * ), NDIMR( * );
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               TWO
      PARAMETER          ( TWO = 2.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I, IL, IR, LLST, MAXN, NCRNT, NLVL;
      REAL               TEMP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, LOG, MAX, REAL
      // ..
      // .. Executable Statements ..
*
      // Find the number of levels on the tree.
*
      MAXN = MAX( 1, N )
      TEMP = LOG( REAL( MAXN ) / REAL( MSUB+1 ) ) / LOG( TWO )
      LVL = INT( TEMP ) + 1
*
      I = N / 2
      INODE( 1 ) = I + 1
      NDIML( 1 ) = I
      NDIMR( 1 ) = N - I - 1
      IL = 0
      IR = 1
      LLST = 1
      DO 20 NLVL = 1, LVL - 1
*
         // Constructing the tree at (NLVL+1)-st level. The number of
         // nodes created on this level is LLST * 2.
*
         DO 10 I = 0, LLST - 1
            IL = IL + 2
            IR = IR + 2
            NCRNT = LLST + I
            NDIML( IL ) = NDIML( NCRNT ) / 2
            NDIMR( IL ) = NDIML( NCRNT ) - NDIML( IL ) - 1
            INODE( IL ) = INODE( NCRNT ) - NDIMR( IL ) - 1
            NDIML( IR ) = NDIMR( NCRNT ) / 2
            NDIMR( IR ) = NDIMR( NCRNT ) - NDIML( IR ) - 1
            INODE( IR ) = INODE( NCRNT ) + NDIML( IR ) + 1
   10    CONTINUE
         LLST = LLST*2
   20 CONTINUE
      ND = LLST*2 - 1
*
      RETURN
*
      // End of SLASDT
*
      END
