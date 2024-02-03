      int     FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK );
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                IHI, ILO, ISPEC, LWORK, N;
      String             NAME*( * ), OPTS*( * );
*
*  ================================================================
*     .. Parameters ..
      int                INMIN, INWIN, INIBL, ISHFTS, IACC22, ICOST;
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14, ISHFTS = 15, IACC22 = 16, ICOST = 17 )
      int                NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP, RCOST;
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14, NIBBLE = 14, KNWSWP = 500, RCOST = 10 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
*     ..
*     .. Local Scalars ..
      int                NH, NS;
      int                I, IC, IZ;
      String             SUBNAM*6;
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC LOG, MAX, MOD, NINT, REAL
*     ..
*     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR. ( ISPEC.EQ.IACC22 ) ) THEN
*
*        ==== Set the number simultaneous shifts ====
*
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 ) NS = 4          IF( NH.GE.60 ) NS = 10          IF( NH.GE.150 ) NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )          IF( NH.GE.590 ) NS = 64          IF( NH.GE.3000 ) NS = 128          IF( NH.GE.6000 ) NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
*
      IF( ISPEC.EQ.INMIN ) THEN
*
*
*        ===== Matrices of order smaller than NMIN get sent
*        .     to xLAHQR, the classic double shift algorithm.
*        .     This must be at least 11. ====
*
         IPARMQ = NMIN
*
      ELSE IF( ISPEC.EQ.INIBL ) THEN
*
*        ==== INIBL: skip a multi-shift qr iteration and
*        .    whenever aggressive early deflation finds
*        .    at least (NIBBLE*(window size)/100) deflations. ====
*
         IPARMQ = NIBBLE
*
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
*
*        ==== NSHFTS: The number of simultaneous shifts =====
*
         IPARMQ = NS
*
      ELSE IF( ISPEC.EQ.INWIN ) THEN
*
*        ==== NW: deflation window size.  ====
*
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
*
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
*
*        ==== IACC22: Whether to accumulate reflections
*        .     before updating the far-from-diagonal elements
*        .     and whether to use 2-by-2 block structure while
*        .     doing it.  A small amount of work could be saved
*        .     by making this choice dependent also upon the
*        .     NH=IHI-ILO+1.
*
*
*        Convert NAME to upper case if the first character is lower case.
*
         IPARMQ = 0
         SUBNAM = NAME
         IC = ICHAR( SUBNAM( 1: 1 ) )
         IZ = ICHAR( 'Z' )
         IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*           ASCII character set
*
            IF( IC.GE.97 .AND. IC.LE.122 ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( IC.GE.97 .AND. IC.LE.122 ) SUBNAM( I: I ) = CHAR( IC-32 )
               END DO
            END IF
*
         ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*           EBCDIC character set
*
            IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. ( IC.GE.145 .AND. IC.LE.153 ) .OR. ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC+64 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. ( IC.GE.145 .AND. IC.LE.153 ) .OR. ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I: I ) = CHAR( IC+64 )
               END DO
            END IF
*
         ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*           Prime machines:  ASCII+128
*
            IF( IC.GE.225 .AND. IC.LE.250 ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( IC.GE.225 .AND. IC.LE.250 ) SUBNAM( I: I ) = CHAR( IC-32 )
               END DO
            END IF
         END IF
*
         IF( SUBNAM( 2:6 ).EQ.'GGHRD' .OR. SUBNAM( 2:6 ).EQ.'GGHD3' ) THEN
            IPARMQ = 1
            IF( NH.GE.K22MIN ) IPARMQ = 2
         ELSE IF ( SUBNAM( 4:6 ).EQ.'EXC' ) THEN
            IF( NH.GE.KACMIN ) IPARMQ = 1             IF( NH.GE.K22MIN ) IPARMQ = 2          ELSE IF ( SUBNAM( 2:6 ).EQ.'HSEQR' .OR. SUBNAM( 2:5 ).EQ.'LAQR' ) THEN             IF( NS.GE.KACMIN ) IPARMQ = 1             IF( NS.GE.K22MIN ) IPARMQ = 2
         END IF
*
      ELSE IF( ISPEC.EQ.ICOST ) THEN
*
*        === Relative cost of near-the-diagonal chase vs
*            BLAS updates ===
*
         IPARMQ = RCOST
      ELSE
*        ===== invalid value of ispec =====
         IPARMQ = -1
*
      END IF
*
*     ==== End of IPARMQ ====
*
      END
