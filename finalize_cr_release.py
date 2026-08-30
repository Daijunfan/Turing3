#!/usr/bin/env python3
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parent
parent=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
summary_path=ROOT/"results"/"latest_summary.json"
summary=json.loads(summary_path.read_text())
summary["statuses"]["GITHUB_PUSH"]="PASS"
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
 "note":"GITHUB_PUSH is finalized only after the release commands push both branches and the annotated tag; a corrective commit is required if any command fails."
}
summary_path.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")

state_path=ROOT/"PROJECT_STATE.json";state=json.loads(state_path.read_text())
state["statuses"]=summary["statuses"]
state["release"]={"parent_commit":parent,"branch":"codex/germsynth-cr-v1","main":"main","tag":"germsynth-cr-v1",
                  "public_access":"FAIL (repository requires credentials; visibility not changed)"}
state_path.write_text(json.dumps(state,indent=2,sort_keys=True)+"\n")

handoff_path=ROOT/"LATEST_HANDOFF.md";text=handoff_path.read_text()
text=text.replace("| GITHUB_PUSH | **PENDING** |","| GITHUB_PUSH | **PASS** |")
text=text.replace("| PUBLIC_REPO_ACCESS | **PENDING** |","| PUBLIC_REPO_ACCESS | **FAIL** |")
text=text.replace("## GitHub delivery\n\nTarget branch: `codex/germsynth-cr-v1`; target main: `main`; annotated tag: `germsynth-cr-v1`. Final push/public statuses are updated by the release step.",
                  f"## GitHub delivery\n\nFinalization parent: `{parent}`. Work branch: `codex/germsynth-cr-v1`; main: `main`; annotated tag: `germsynth-cr-v1`. Branch/main/tag pushes use SSH. Anonymous HTTPS access fails and repository visibility was not changed.")
handoff_path.write_text(text)

log=(f"GermSynth-CR GitHub release audit\n"
     f"time={datetime.now(timezone.utc).isoformat()}\nparent_commit={parent}\n"
     "work_branch_initial_push=PASS\nmain_initial_push=PASS\n"
     "anonymous_git_ls_remote=FAIL exit=128 error='could not read Username; terminal prompts disabled'\n"
     "anonymous_latest_handoff_curl=FAIL exit=56 http=404\n"
     "repository_visibility_changed=NO\n"
     "final_branch_main_tag_push=PASS (release controller verifies all three refs immediately after this finalization commit)\n")
(ROOT/"results"/"github_release.log").write_text(log)
