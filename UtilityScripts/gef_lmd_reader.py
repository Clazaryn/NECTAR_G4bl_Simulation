# gef_lmd_reader.py  (Python 3.5 compatible)
import os
import sys
import time


""" USING THIS CLASS YOU WILL GET WHATEVER YOU WANT FROM THE .lmd file of GEF 
    assuming your GEF folder is next to the one of NECTAR and that you already generated
    the .lmd file with GEFbashscript.sh

    You can test and see what beautiful plot you can do while running the GEFexample_1.py launch from the Nectar folder 

"""


class GEFEvent(object):
    """
    Simple container for the 27 fields read in the .lmd event line.
    """

    def __init__(self,
                 FissionMode,
                 Zsad,
                 Asad,
                 Qvalue,
                 Z_light,
                 Z_heavy,
                 A_light_pre,
                 A_heavy_pre,
                 A_light_post,
                 A_heavy_post,
                 I1pre,
                 I2pre,
                 I1gs,
                 I2gs,
                 T_light_pre,
                 cos_theta_light,
                 phi_light_deg,
                 T_heavy_pre,
                 cos_theta_heavy,
                 phi_heavy_deg,
                 Eexc1,
                 Eexc2,
                 n1,
                 n2,
                 T_light_post,
                 T_heavy_post,
                 E_at_fission):

        self.FissionMode = FissionMode
        self.Zsad = Zsad
        self.Asad = Asad
        self.Qvalue = Qvalue
        self.Z_light = Z_light
        self.Z_heavy = Z_heavy
        self.A_light_pre = A_light_pre
        self.A_heavy_pre = A_heavy_pre
        self.A_light_post = A_light_post
        self.A_heavy_post = A_heavy_post
        self.I1pre = I1pre
        self.I2pre = I2pre
        self.I1gs = I1gs
        self.I2gs = I2gs
        self.T_light_pre = T_light_pre
        self.cos_theta_light = cos_theta_light
        self.phi_light_deg = phi_light_deg
        self.T_heavy_pre = T_heavy_pre
        self.cos_theta_heavy = cos_theta_heavy
        self.phi_heavy_deg = phi_heavy_deg
        self.Eexc1 = Eexc1
        self.Eexc2 = Eexc2
        self.n1 = n1
        self.n2 = n2
        self.T_light_post = T_light_post
        self.T_heavy_post = T_heavy_post
        self.E_at_fission = E_at_fission


class GEFLmdReader(object):
    """
    Reads a GEF .lmd (ASCII) file without ROOT.

    New signature you asked:
        GEFLmdReader(Z, A, Estar, factor, event_number, base_dir=".")

    - event_number:
        If not None, only builds the index up to this many events (faster for huge files).
        If None or <=0, indexes the full file.

    Random access uses a byte-offset index, so once indexed:
        get_event(i) is O(1) in file size.
    """

    def __init__(self,
                 Z,
                 A,
                 Estar,
                 factor,
                 event_number=None,
                 base_dir=".",
                 filename_template="../GEF/out/GEFResults_Z{Z}_A{A}_E{E}_factor_{factor}.lmd",
                 one_based_events=True,
                 build_index=True,
                 show_progress=True):

        self.Z = int(Z)
        self.A = int(A)
        self.Estar = float(Estar)
        self.factor = int(factor)

        self.base_dir = base_dir
        self.filename_template = filename_template
        self.one_based_events = bool(one_based_events)
        self.show_progress = bool(show_progress)

        # Limit for indexing (None => full file)
        self.max_events = None
        if event_number is not None:
            self.max_events = int(event_number)
            if self.max_events <= 0:
                self.max_events = None

        self.path = self._resolve_file()

        self._offsets = []
        self._fh = None  # optional persistent handle (use open()/close())

        if build_index:
            self._build_index(max_events=self.max_events)

    # ---------- file selection ----------
    def _resolve_file(self):
        # match filenames like ..._E2.0_factor_1.lmd
        e_str = "{0:.1f}".format(self.Estar)

        fname = self.filename_template.format(
            Z=self.Z,
            A=self.A,
            E=e_str,
            factor=self.factor
        )

        p = os.path.abspath(os.path.join(self.base_dir, fname))
        if not os.path.exists(p):
            raise IOError(
                "GEFLmdReader ERROR: file not found:\n{0}\n"
                "Check Z, A, E*, factor or base_dir/filename_template.".format(p)
            )
        return p

    # ---------- open/close for fast repeated random access ----------
    def open(self):
        if self._fh is None:
            try:
                self._fh = open(self.path, "r")
            except IOError as e:
                raise IOError(
                    "GEFLmdReader ERROR: cannot open file '{0}'.\n"
                    "Original error: {1}".format(self.path, e)
                )
        return self

    def close(self):
        if self._fh is not None:
            try:
                self._fh.close()
            finally:
                self._fh = None

    # ---------- indexing helpers ----------
    def _is_data_line(self, line):
        s = line.lstrip()
        if not s:
            return False
        c0 = s[0]
        if c0.isdigit():
            return True
        if c0 == "-" and len(s) > 1 and s[1].isdigit():
            return True
        return False

    def _try_parse_line(self, line):
        # C++ reader reads exactly 27 fields; ignore extra trailing fields
        parts = line.split()
        if len(parts) < 27:
            return None
        try:
            return GEFEvent(
                int(parts[0]),     # FissionMode
                int(parts[1]),     # Zsad
                int(parts[2]),     # Asad
                float(parts[3]),   # Qvalue
                int(parts[4]),     # Z_light
                int(parts[5]),     # Z_heavy
                int(parts[6]),     # A_light_pre
                int(parts[7]),     # A_heavy_pre
                int(parts[8]),     # A_light_post
                int(parts[9]),     # A_heavy_post
                float(parts[10]),  # I1pre
                float(parts[11]),  # I2pre
                float(parts[12]),  # I1gs
                float(parts[13]),  # I2gs
                float(parts[14]),  # T_light_pre
                float(parts[15]),  # cos_theta_light
                float(parts[16]),  # phi_light_deg
                float(parts[17]),  # T_heavy_pre
                float(parts[18]),  # cos_theta_heavy
                float(parts[19]),  # phi_heavy_deg
                float(parts[20]),  # Eexc1
                float(parts[21]),  # Eexc2
                int(parts[22]),    # n1
                int(parts[23]),    # n2
                float(parts[24]),  # T_light_post
                float(parts[25]),  # T_heavy_post
                float(parts[26])   # E_at_fission
            )
        except Exception:
            return None

    def _build_index(self, max_events=None):
        """
        Build offset index.
        If max_events is not None: stop after indexing max_events events.
        Shows a progress bar based on bytes read.
        """
        self._offsets = []

        total_bytes = None
        if self.show_progress:
            try:
                total_bytes = os.path.getsize(self.path)
            except Exception:
                total_bytes = None

        t0 = time.time()
        last_print = 0.0

        f = open(self.path, "r")
        try:
            while True:
                pos = f.tell()
                line = f.readline()
                if not line:
                    break

                if self._is_data_line(line):
                    ev = self._try_parse_line(line)
                    if ev is not None:
                        self._offsets.append(pos)

                        if (max_events is not None) and (len(self._offsets) >= max_events):
                            break

                if self.show_progress:
                    now = time.time()
                    if (now - last_print) > 0.25:
                        last_print = now

                        done = f.tell()
                        elapsed = now - t0
                        speed = (done / 1024.0 / 1024.0) / elapsed if elapsed > 0 else 0.0  # MB/s

                        # If we are stopping early, it's nicer to show event-based progress:
                        if max_events is not None and max_events > 0:
                            frac = float(len(self._offsets)) / float(max_events)
                            if frac > 1.0:
                                frac = 1.0
                            bar_len = 30
                            filled = int(bar_len * frac)
                            bar = "=" * filled + "-" * (bar_len - filled)
                            sys.stdout.write(
                                "\rIndexing [{0}] {1:6.2f}%  ({2:.1f} MB/s)  events={3}/{4}".format(
                                    bar, 100.0 * frac, speed, len(self._offsets), max_events
                                )
                            )
                        else:
                            # Full file: byte-based progress
                            if total_bytes and total_bytes > 0:
                                frac = float(done) / float(total_bytes)
                                if frac > 1.0:
                                    frac = 1.0
                                bar_len = 30
                                filled = int(bar_len * frac)
                                bar = "=" * filled + "-" * (bar_len - filled)
                                sys.stdout.write(
                                    "\rIndexing [{0}] {1:6.2f}%  ({2:.1f} MB/s)  events={3}".format(
                                        bar, 100.0 * frac, speed, len(self._offsets)
                                    )
                                )
                            else:
                                sys.stdout.write(
                                    "\rIndexing... ({0:.1f} MB/s)  events={1}".format(speed, len(self._offsets))
                                )

                        sys.stdout.flush()

        finally:
            f.close()

        if self.show_progress:
            sys.stdout.write("\nDone indexing. Found {0} events.\n".format(len(self._offsets)))
            sys.stdout.flush()

    # ---------- event access ----------
    def _event_index(self, event_number):
        n = int(event_number)
        if self.one_based_events:
            idx = n - 1
        else:
            idx = n

        if idx < 0 or idx >= len(self._offsets):
            raise IndexError("Event {0} out of range. Valid: 1..{1}".format(
                event_number, len(self._offsets)
            ))
        return idx

    def n_events(self):
        return len(self._offsets)

    def get_event(self, event_number):
        idx = self._event_index(event_number)

        if self._fh is None:
            f = open(self.path, "r")
            try:
                f.seek(self._offsets[idx])
                line = f.readline()
            finally:
                f.close()
        else:
            self._fh.seek(self._offsets[idx])
            line = self._fh.readline()

        ev = self._try_parse_line(line)
        if ev is None:
            raise RuntimeError("Failed to parse event line for event {0}".format(event_number))
        return ev

    # ---------- convenience getters ----------
    def GetEexc1(self, event_number):
        return self.get_event(event_number).Eexc1

    def GetEexc2(self, event_number):
        return self.get_event(event_number).Eexc2

    def GetZ_pair(self, event_number):
        ev = self.get_event(event_number)
        return ev.Z_light, ev.Z_heavy

    def GetApre_pair(self, event_number):
        ev = self.get_event(event_number)
        return ev.A_light_pre, ev.A_heavy_pre

    def GetApost_pair(self, event_number):
        ev = self.get_event(event_number)
        return ev.A_light_post, ev.A_heavy_post

    def GetKEpre_pair(self, event_number):
        ev = self.get_event(event_number)
        return ev.T_light_pre, ev.T_heavy_pre

    def GetKEpost_pair(self, event_number):
        ev = self.get_event(event_number)
        return ev.T_light_post, ev.T_heavy_post

    def GetNeutrons_pair(self, event_number):
        ev = self.get_event(event_number)
        return ev.n1, ev.n2
    
    def GetCosTheta_pair(self, event_number):
        ev = self.get_event(event_number)
        return ev.cos_theta_light, ev.cos_theta_heavy

    def GetPhi_pair(self, event_number):
        ev = self.get_event(event_number)
        return ev.phi_light_deg, ev.phi_heavy_deg
