# mose_completion.sh -- Bash completion for the mose launcher. Add to ~/.bashrc:
#   source /path/to/MOSE-GPU/scripts/mose_completion.sh
# Tab then completes the commands, the options of each command, the case folders, and the SLURM partitions after -p.

_mose_completions()
{
  local cur prev cmd i
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"
  cmd=""
  for ((i = 1; i < COMP_CWORD; i++)); do
    case ${COMP_WORDS[i]} in
      check|build|solve|submit|list|example) cmd=${COMP_WORDS[i]}; break ;;
    esac
  done

  case ${prev} in
    -p|--partition)
      COMPREPLY=( $(compgen -W "$(sinfo -h -o '%P %G' 2>/dev/null | awk '$2 ~ /gpu/ {sub(/\*$/, "", $1); print $1}' | sort -u)" -- "${cur}") )
      return 0 ;;
    --build-partition)
      COMPREPLY=( $(compgen -W "$(sinfo -h -o %P 2>/dev/null | tr -d '*' | sort -u)" -- "${cur}") )
      return 0 ;;
    --nodelist)
      COMPREPLY=( $(compgen -W "$(sinfo -h -o %n 2>/dev/null | sort -u)" -- "${cur}") )
      return 0 ;;
    --src|--builds)
      compopt -o filenames 2>/dev/null
      COMPREPLY=( $(compgen -d -- "${cur}") )
      return 0 ;;
    -n|--ntasks|--gpus|-t|--time|-J|--name|--mem|--cpus-per-task)
      return 0 ;;
  esac

  case ${cmd} in
    "")
      COMPREPLY=( $(compgen -W "check build solve list example --src --builds --build-partition -h" -- "${cur}") ) ;;
    solve|submit|check|build)
      if [[ ${cur} == -* ]]; then
        if [[ ${cmd} == solve || ${cmd} == submit ]]; then
          COMPREPLY=( $(compgen -W "-n --gpus -p -t -J --mem --cpus-per-task --nodelist --dry-run -h" -- "${cur}") )
        else
          COMPREPLY=( $(compgen -W "--dry-run -h" -- "${cur}") )
        fi
      else
        compopt -o filenames 2>/dev/null
        COMPREPLY=( $(compgen -d -- "${cur}") )
      fi ;;
  esac
  return 0
}

complete -F _mose_completions mose
