import straight.plugin


def load_plugins():
    '''Load all available types of plugins'''
    plugins = {
        'cluster': [],
        'generic': [],
        'genome_wide': [],
        'input': [],
        'output': [],
        'specific': []
    }

    for key in plugins.keys():
        plugins[key] = list(straight.plugin.load('antismash.{}'.format(key)))
        plugins[key].sort(key=lambda x: x.priority)

    return plugins


def check_prereqs(plugins, options):
    '''Check if all prerequisites are met'''
    missing_prereqs = []

    for _, plugin_list in plugins.iteritems():
        for plugin in plugin_list:
            missing_prereqs.extend(plugin.check_prereqs(options))

    # remove duplicates
    missing_prereqs = list(set(missing_prereqs))
    missing_prereqs.sort()
    return missing_prereqs


def get_versions(plugins, options):
    versions = []
    for _, plugin_list in plugins.iteritems():
        for plugin in plugin_list:
            versions.extend(plugin.get_versions(options))

    # remove duplicate entries
    versions = list(set(versions))
    versions.sort()
    return versions
