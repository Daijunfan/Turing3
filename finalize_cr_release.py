#!/usr/bin/env python3
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parent
parent=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
summary_path=ROOT/"results"/"latest_summary.json"
summary=json.loads(summary_path.read_text())
summary["statuses"]["GITHUB_PUSH"]="FAIL"
summary["statuses"]["PUBLIC_REPO_ACCESS"]="FAIL"
summary["github_release"]={
 "finalization_parent_commit":parent,
 "work_branch":"codex/germsynth-cr-v1",
 "main":"main",
 "annotated_tag":"germsynth-cr-v1",
 "push_transport":"SSH",
 "anonymous_git_ls_remote":{"status":"FAIL","exit":128,"error":"could not read Username; terminal prompts disabled"},
 "anonymous_raw_handoff":{"status":"FAIL","curl_exit":56,"http":"404"},
 "visibility_action":"UNCHANGED",
 "note":"SSH branch/main/tag pushes passed, but the project gate requires anonymous handoff access too; therefore GITHUB_PUSH and PUBLIC_REPO_ACCESS are both FAIL."
}
summary_path.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")

state_path=ROOT/"PROJECT_STATE.json";state=json.loads(state_path.read_text())
state["base_commit"]=parent
state["statuses"]=summary["statuses"]
state["release"]={"parent_commit":parent,"branch":"codex/germsynth-cr-v1","main":"main","tag":"germsynth-cr-v1",
                  "public_access":"FAIL (repository requires credentials; visibility not changed)"}
state_path.write_text(json.dumps(state,indent=2,sort_keys=True)+"\n")

handoff_path=ROOT/"LATEST_HANDOFF.md";text=handoff_path.read_text()
text=re.sub(r"Generated from commit: `[^`]+` \(the delivery tag is authoritative because a commit cannot contain its own SHA\)\.",
            f"Release-state parent commit: `{parent}`. The checked-out `main` SHA is authoritative because a commit cannot contain its own SHA.", text)
text=text.replace("| GITHUB_PUSH | **PENDING** |","| GITHUB_PUSH | **FAIL** |")
text=text.replace("| GITHUB_PUSH | **PASS** |","| GITHUB_PUSH | **FAIL** |")
text=text.replace("| PUBLIC_REPO_ACCESS | **PENDING** |","| PUBLIC_REPO_ACCESS | **FAIL** |")
start=text.index("## GitHub delivery")
text=text[:start] + ("## GitHub delivery\n\n"
    f"Release-state parent: `{parent}`. Work branch and main contain this corrective state after synchronization. "
    "Annotated tag `germsynth-cr-v1` points to the preceding release commit `9e3ebe640f373ba91cd60eaef4369daf5cd179f8` and is not force-moved. "
    "SSH branch/main/tag pushes passed, but the combined GITHUB gate is FAIL because anonymous HTTPS access fails; repository visibility was not changed.\n")
text=text.replace("## GitHub delivery\n\nTarget branch: `codex/germsynth-cr-v1`; target main: `main`; annotated tag: `germsynth-cr-v1`. Final push/public statuses are updated by the release step.",
                  f"## GitHub delivery\n\nFinalization parent: `{parent}`. Work branch and main are updated through the corrective commit. Annotated tag `germsynth-cr-v1` remains on the preceding release commit because force-moving tags is forbidden. SSH pushes passed, but the combined GITHUB gate is FAIL because anonymous HTTPS access fails; repository visibility was not changed.")
handoff_path.write_text(text)

log=(f"GermSynth-CR GitHub release audit\n"
     f"time={datetime.now(timezone.utc).isoformat()}\nparent_commit={parent}\n"
     "work_branch_initial_push=PASS\nmain_initial_push=PASS\n"
     "anonymous_git_ls_remote=FAIL exit=128 error='could not read Username; terminal prompts disabled'\n"
     "anonymous_latest_handoff_curl=FAIL exit=56 http=404\n"
     "repository_visibility_changed=NO\n"
     "ssh_branch_main_tag_push=PASS\ncombined_GITHUB_PUSH_gate=FAIL (anonymous access is required and failed)\n")
(ROOT/"results"/"github_release.log").write_text(log)
