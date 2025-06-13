def format_seconds(seconds):
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    
    if hours > 0:
        return f"{hours:02d}:{minutes:02d}:{secs:02d} hours"
    elif minutes > 0:
        return f"{minutes:02d}:{secs:02d} minutes"
    else:
        return f"{secs:02d} seconds"