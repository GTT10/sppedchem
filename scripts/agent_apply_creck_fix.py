from pathlib import Path
import re


def apply() -> None:
    source_path = Path("src/chemkin_module.f90")
    source = source_path.read_text()

    if "Continuing with II pinned at IDIM" not in source:
        pattern = re.compile(
            r"(?m)^            ELSE\n"
            r"               WRITE \(LOUT, 1070\)\n"
            r"               KERR = \.TRUE\.\n"
            r"            ENDIF$"
        )
        replacement = """            ELSE
!              Fail immediately. Continuing with II pinned at IDIM makes
!              auxiliary records from later reactions attach to reaction
!              IDIM and can replace this primary error with a false PLOG
!              ordering diagnostic.
               WRITE (*,1070) IDIM
               WRITE (LOUT,1070) IDIM
               CLOSE (LIN)
               CLOSE (LOUT)
               ERROR STOP 1
            ENDIF"""
        source, count = pattern.subn(replacement, source, count=1)
        if count != 1:
            raise RuntimeError(f"reaction overflow block matches: {count}")

    if "Preserve the primary parser diagnostic" not in source:
        search_from = source.index(" 5000 CONTINUE")
        start = source.index("      CLOSE (LIN)\n", search_from)
        end = source.index("      DO 1150 K = 1, KK\n", start)
        replacement = """      CLOSE (LIN)
!
!     Preserve the primary parser diagnostic. Do not run PLOG validation
!     and do not create a linking file when an earlier parse check failed.
      IF (KERR) THEN
         WRITE (LOUT, '(//A)')&
         ' ERROR...MECHANISM PARSING FAILED; NO LINKING FILE WAS WRITTEN'
         CLOSE (LOUT)
         ERROR STOP 1
      ENDIF
!
!     Validate and pack PLOG data before opening cklink. plog_finalize uses
!     ERROR STOP for malformed data, so no partial link exists on failure.
      CALL plog_finalize(KERR, LOUT)
      IF (KERR) THEN
         WRITE (LOUT, '(//A)')&
         ' ERROR...PLOG VALIDATION FAILED; NO LINKING FILE WAS WRITTEN'
         CLOSE (LOUT)
         ERROR STOP 1
      ENDIF
!
!     Create the link only after all parse-time validation succeeds.
      OPEN (LINC, FORM='UNFORMATTED', STATUS='REPLACE',&
                  FILE=trim(mechdir)//"cklink")
      REWIND LINC
!     v2 leading record: magic + integer schema version. Positively
!     identifies the format for the reader (SCcklink) — supersedes the
!     fragile VERS-string check.
      WRITE (LINC) CK_MAGIC, CK_SCHEMA
      WRITE (LINC) VERS, PREC, KERR

!
"""
        source = source[:start] + replacement + source[end:]

    if "reaction count exceeds CKINTP capacity IDIM=" not in source:
        pattern = re.compile(
            r"(?m)^ 1070 FORMAT \(6X,'Error\.\.\.more than IDIM reactions\.\.\.'\)$"
        )
        replacement = (
            " 1070 FORMAT (6X,'ERROR...reaction count exceeds CKINTP capacity IDIM=',&\n"
            "                   I0,'. Aborting before auxiliary/PLOG parsing.')"
        )
        source, count = pattern.subn(replacement, source, count=1)
        if count != 1:
            raise RuntimeError(f"format 1070 matches: {count}")

    source_path.write_text(source)


apply()
