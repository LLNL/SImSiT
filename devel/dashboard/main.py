import os
import dateutil
from datetime import date
from random import randint
import numpy as np

from bokeh.plotting import figure, output_file, show, curdoc
from bokeh.layouts import grid, gridplot, column, row
from bokeh.models.widgets import Panel, Tabs, DataTable, DateFormatter, TableColumn
from bokeh.io import output_file, show, curdoc
from bokeh.models import ColumnDataSource, PreText, Select, CustomJS, DatetimeTickFormatter, Div, LabelSet
from bokeh.models.tools import HoverTool
from bokeh.core.properties import FontSize
from bokeh.core.enums import FontStyle
from bokeh.transform import cumsum

#output to static HTML file
#output_file("xfiles.html")

k_figtools = "pan,reset,save"

def parseDirectories():
    """Goes through directories and places the data into dictionaries based on branch and competitor
    """
    path_to_branch_names = 'submissions/'
    branches = os.listdir(path=path_to_branch_names)
    try:
        branches.remove('.DS_Store')
    except ValueError:
        pass

    parsed_data = {}
    for branch in branches:
        branch_path = os.path.join(path_to_branch_names, branch)
        parsed_data[branch] = {}
        combined_data = []
        for file in os.listdir(branch_path):
            if file.endswith('.txt'):
                file_name = os.path.join(branch_path, file)
                with open(file_name) as fd:
                    for line in fd:
                        line = line.split()
                        competitor = line[0].lower()
                        parsed_data[branch][competitor] = parsed_data[branch].get(competitor, [])
                        parsed_data[branch][competitor].append(line[2:])
                        combined_data.append(line)
        parsed_data[branch]['summary'] = combined_data
    return parsed_data

def createSummaryData(data):
    """Package data for display on the 'Summary' tab

    Parameters
    ----------
    data : dict
        The output of the `parseDirectories` function with data loaded from the score tables
    """

    branches = ['iod', 'detect_sidereal', 'detect_target', 'calibrate_sidereal']
    data['summary'] = []
    for branch in branches:
        comp = []
        date = []
        cur_max_competitor = ''
        score = []
        cur_max_score = None
        for entry in sorted(data[branch]['summary'], key=lambda x: x[2]):
            date.append(dateutil.parser.parse(entry[2]))
            if branch == 'iod':
                if (cur_max_score is None) or (float(entry[10]) > cur_max_score):
                    cur_max_competitor = makePrettyName(entry[0])
                    cur_max_score = float(entry[10])
            else:
                if (cur_max_score is None) or (float(entry[4])/float(entry[3]) > cur_max_score):
                    cur_max_competitor = makePrettyName(entry[0])
                    cur_max_score = float(entry[4])/float(entry[3])
            comp.append(cur_max_competitor)
            score.append(cur_max_score)
        data['summary'] += [[branch, comp, date, score]]

    for branch in branches:
        data[branch]['summary_table'] = {'date':[], 'score':[], 'competitors': []}
        for competitor, values in data[branch].items():
            if competitor != 'summary_table' and competitor != 'summary':
                max_score = None
                for value in values:
                    if branch == 'iod':
                        if max_score is None or float(value[8]) > float(max_score[8]):
                            max_score = value
                    else:
                        if max_score is None or float(value[2])/float(value[1]) > float(max_score[2])/float(max_score[1]):
                            max_score = value
                if branch == 'iod':
                    data[branch]['summary_table']['date'].append(dateutil.parser.parse(max_score[0]))
                    data[branch]['summary_table']['score'].append(float(max_score[8]))
                    data[branch]['summary_table']['competitors'].append(makePrettyName(competitor))
                else:
                    data[branch]['summary_table']['date'].append(dateutil.parser.parse(max_score[0]))
                    data[branch]['summary_table']['score'].append(float(max_score[2])/float(max_score[1]))
                    data[branch]['summary_table']['competitors'].append(makePrettyName(competitor))

    return data

def makeKeyName(name):
    """Translate a display name into an easily parsed key name for dictionary indexing
    """
    return '_'.join(name.split()).lower()

def makePrettyName(name):
    """Translate a key value into a nicely formatted display name
    """
    if name == 'iod':
        return "IOD"
    spaced_name = name.split('_')
    return ' '.join(list(map(lambda x: x[0].upper() + x[1:], spaced_name)))

def onTabChange(attr, old, new):
    """Execute these commands when a new top-level tab is selected
    """
    current_tab = makeKeyName(tabs.tabs[new].title)
    if current_tab == "summary":
        current_competitor = None
    else:
        select_options = list(data[current_tab].keys())
        try:
            select_options.remove('summary')
            select_options.remove('summary_table')
        except ValueError:
            pass
        select_bar.options = list(map(makePrettyName, select_options))
        current_competitor = select_bar.options[0]
    select_bar.value = current_competitor

    redrawPlot(current_tab, current_competitor)

def onSummarySelectChange(attr, old, new):
    data_sources['summary_table'].data = data[makeKeyName(new)]['summary_table']

def onSelectChange(attr, old, new):
    """Execute these commands when the pull-down menu on branch tabs is adjusted
    """
    current_tab = tabs.tabs[tabs.active].title
    current_competitor = select_bar.value

    redrawPlot(current_tab, current_competitor)

def redrawPlot(tab, competitor):
    """Update the data that is rendered in each plot

    This method is called each time a tab or competitor selection is changed to refresh the data 
    for plotting.
    """
    tab = makeKeyName(tab)
    for row_idx, row in enumerate(display_shape[tab]):
        for idx, (title, (x_axis_label, y_axis_label), plot, data_idx) in enumerate(row):
            if '-' in data_idx:
                data_idx = data_idx.split('-')
            else:
                data_idx = int(data_idx) + 1

            colors = []
            x = []
            y = []
            competitor = makeKeyName(competitor)
            for submission in data[tab][competitor]:
                if plot == "xy":
                    x.append(dateutil.parser.parse(submission[0]))
                    y.append(float(submission[data_idx]))
                    data_sources[tab][row_idx][idx].data = {'x': x, 'y': y}
                elif plot == "xy_history":
                    cur_y = float(submission[data_idx])*100
                    x.append(dateutil.parser.parse(submission[0]))
                    y.append(cur_y)
                    if cur_y >= 70:
                        colors.append('green')
                    elif cur_y < 70 and cur_y >= 40:
                        colors.append('gold')
                    else:
                        colors.append('red')
                    data_sources[tab][row_idx][idx].data = {'x': x, 'y': y, 'colors': colors}
                elif plot == "xy_history_div":
                    det_scores = list(map(float, submission[int(data_idx[0])+1:int(data_idx[1])+2]))
                    cur_y = (det_scores[1]/det_scores[0])*100
                    x.append(dateutil.parser.parse(submission[0]))
                    y.append(cur_y)
                    if cur_y >= 70:
                        colors.append('green')
                    elif cur_y < 70 and cur_y >= 40:
                        colors.append('gold')
                    else:
                        colors.append('red')
                    data_sources[tab][row_idx][idx].data = {'x': x, 'y': y, 'colors': colors}
                elif plot == "xy_history_div_fpr":
                    det_scores = list(map(float, submission[int(data_idx[0])+1:int(data_idx[1])+2]))
                    cur_y = (1-(det_scores[1]/det_scores[0]))*100
                    x.append(dateutil.parser.parse(submission[0]))
                    y.append(cur_y)
                    if cur_y >= 60:
                        colors.append('red')
                    elif cur_y < 60 and cur_y >= 30:
                        colors.append('gold')
                    else:
                        colors.append('green')
                    data_sources[tab][row_idx][idx].data = {'x': x, 'y': y, 'colors': colors}
                elif plot == "bar":
                    colors = []
                    det_scores = list(map(float, submission[int(data_idx[0])+1:int(data_idx[1])+2]))
                    y = [det_scores[1]/det_scores[0], 1-(det_scores[3]/det_scores[2])]
                    y = list(map(lambda x: x * 100.0, y))
                    if y[0] >= 70:
                        colors.append('green')
                    elif y[0] < 70 and y[0] >= 40:
                        colors.append('gold')
                    else:
                        colors.append('red')

                    if y[1] >= 60:
                        colors.append('red')
                    elif y[1] < 60 and y[1] >= 30:
                        colors.append('gold')
                    else:
                        colors.append('green')

                    data_sources[tab][row_idx][idx].data = {'x': [.5, 1.5], 'y': y, 'colors': colors}
                elif plot == "donut":
                    y = float(submission[data_idx])
                    if y >= .70:
                        color = 'green'
                    elif y < .70 and y >= .40:
                        color = 'gold'
                    else:
                        color = 'red'

                    x = [0, 1]
                    y_cur = [np.pi-(y*np.pi), y*np.pi]

                    if y == 1:
                        y_label = [f'{int(y*100):d}%', '']
                    else:
                        y_label = [f'{y*100:.1f}%', '']

                    data_sources[tab][row_idx][idx].data = {'x': x, 'y': y_cur, 'y_label': y_label, 'colors': ['grey', color]} 

def summary_time_series():
    """Plot submission scores over time for each of the branches

    To be displayed on the Summary landing page
    """
    line_colors = ['blue', 'green', 'gold', 'red']
    plot = figure(plot_width=1000, y_range=(0, 1), tools="pan,wheel_zoom,box_zoom,reset", toolbar_location="right", x_axis_type='datetime', 
        title="X-Files submission score history")
    plot.title.text_font_size = "16px"
    plot.xaxis.axis_label = "Date"
    plot.yaxis.axis_label = "Score"
    plot.xaxis.axis_label_text_font_size = "15px"
    plot.yaxis.axis_label_text_font_size = "15px"
    tooltips = [('Date', '@x{%m/%d/%y}'), ('Score', '@y{0.000}'), ('Competitor', '@competitor')]
    hover = HoverTool(tooltips=tooltips, formatters={'@x': 'datetime'})
    plot.add_tools(hover)
    for idx, (branch, comp, date, score) in enumerate(data["summary"]):
        source = ColumnDataSource(data=dict(x=date, y=score, competitor=comp))
        plot.line('x', 'y', source=source, color=line_colors[idx], line_width=3, line_alpha=0.4, legend_label=makePrettyName(branch))
        plot.circle('x', 'y', source=source, size=8, color=line_colors[idx], legend_label=makePrettyName(branch))
    plot.legend.location = "top_left"
    plot.legend.click_policy="hide"
    return [plot]

def leaderboard():
    start_branch = 'iod'
    leaderboard_data = dict(
        date=data[start_branch]['summary_table']['date'],
        score=data[start_branch]['summary_table']['score'],
        competitors=data[start_branch]['summary_table']['competitors']
    )
    
    source = ColumnDataSource(leaderboard_data)

    data_sources['summary_table'] = source

    columns = [
            TableColumn(field="date", title="Date", formatter=DateFormatter()),
            TableColumn(field="score", title="Score"),
            TableColumn(field="competitors", title="Competitor"),
        ]
    data_table = DataTable(source=source, columns=columns, height=280, width=1000)
    return [data_table]

# =============================================================================
# Render the dashboard
# =============================================================================
data = parseDirectories()
data = createSummaryData(data)

# The dashboard starts with the 'Summary' tab selected.
# This summary page is our overview of the rankings across all branches and submission history.
current_tab = "summary"

# On the pages for individual branches, there is a pull-down menu to select which 'Competitor' 
# results to display. There are many performance metrics per branch, which are displayed together 
# on one branch-specific page for a selected competitor.
select_bar = Select(title="Competitor:")
select_bar.on_change('value', onSelectChange)

# The 'tabs' dict holds the summary and branch-specific page specification - one tab for each
tabs = {}
data_sources = {}

# Hard-coded object to represent layout order on webpage.
# List of lists describing the physical layout of each plot. Inisde each list is a tuple with the 
# form ('Title', ('x axis label', 'y axis label'), 'plot type', 'integer index of data to use,
# or a start_index-end_index range') that specifies formatting for individual plot frames.
display_shape = {
            'summary': [],
            'iod': [[('RMSE on x', ('Submission Date', 'RMSE'), 'xy', '0'),
                            ('RMSE on y', ('Submission Date', 'RMSE'), 'xy', '1'),
                            ('RMSE on z', ('Submission Date', 'RMSE'), 'xy', '2')],
                           [('RMSE on vx', ('Submission Date', 'RMSE'), 'xy', '3'),
                            ('RMSE on vy', ('Submission Date', 'RMSE'), 'xy', '4'),
                            ('RMSE on vz', ('Submission Date', 'RMSE'), 'xy',' 5')],
                           [('% Predictions Within FOV After 5 mins', ('', ''), 'donut', '10'),
                            ('% Predictions Within FOV After 15 mins', ('', ''), 'donut', '11'),
                            ('% Predictions Within FOV After 30 mins', ('', ''), 'donut', '12'), 
                            ('% Predictions Within FOV After 60 mins', ('', ''), 'donut', '13')],
                           [('Score History: 5 min', ('Submission Date', 'Score'), 'xy_history', '10'),
                            ('Score History: 15 min', ('Submission Date', 'Score'), 'xy_history', '11'),
                            ('Score History: 30 min', ('Submission Date', 'Score'), 'xy_history', '12'),
                            ('Score History: 60 min', ('Submission Date', 'Score'), 'xy_history', '13')],
                           [('% Pred. in FOV - 5 min (Standardized)', ('', ''), 'donut', '6'),
                            ('% Pred. in FOV - 15 min (Standardized)', ('', ''), 'donut', '7'),
                            ('% Pred. in FOV - 30 min (Standardized)', ('', ''), 'donut', '8'), 
                            ('% Pred. in FOV - 60 min (Standardized)', ('', ''), 'donut', '9')],
                           [('Score History: 5 min - (Standardized)', ('Submission Date', 'Score'), 'xy_history', '6'),
                            ('Score History: 15 min - (Standardized)', ('Submission Date', 'Score'), 'xy_history', '7'),
                            ('Score History: 30 min - (Standardized)', ('Submission Date', 'Score'), 'xy_history', '8'),
                            ('Score History: 60 min - (Standardized)', ('Submission Date', 'Score'), 'xy_history', '9')],
                           [('State Estimate Accuracy', ('Submission Date', 'CVM'), 'xy', '18')]],
            'detect_sidereal': [[('', ('', 'Percent'), 'bar', '0-3'),
                                 ('Completeness History', ('Submission Date', 'Completeness'), 'xy_history_div', '0-1'),
                                 ('False Positive Rate History', ('Submission Date', 'False Positive Rate'), 'xy_history_div_fpr', '2-3')],
                                [('RMSE on position', ('Submission Date', 'RMSE'), 'xy', '4'),
                                 ('RMSE on flux', ('Submission Date', 'RMSE'), 'xy', '5'),
                                 ('RMSE on magnitude', ('Submission Date', 'RMSE'), 'xy', '6')]],
            'detect_target': [[('', ('', 'Percent'), 'bar', '0-3'),
                               ('Completeness History', ('Submission Date', 'Completeness'), 'xy_history_div', '0-1'),
                               ('False Positive Rate History', ('Submission Date', 'False Positive Rate'), 'xy_history_div_fpr', '2-3')],
                              [('RMSE on position', ('Submission Date', 'RMSE'), 'xy', '4'),
                               ('RMSE on flux', ('Submission Date', 'RMSE'), 'xy', '5'),
                               ('RMSE on magnitude', ('Submission Date', 'RMSE'), 'xy', '6')]],
            'calibrate_sidereal': [[('', ('', 'Percent'), 'bar', '0-3'),
                                 ('Completeness History', ('Submission Date', 'Completeness'), 'xy_history_div', '0-1'),
                                 ('False Positive Rate History', ('Submission Date', 'False Positive Rate'), 'xy_history_div_fpr', '2-3')],
                                [('RMSE on position', ('Submission Date', 'RMSE'), 'xy', '4'),
                                 ('RMSE on magnitude', ('Submission Date', 'RMSE'), 'xy', '5')]]
        }

# Create all plots, based on the 'display_shape' dict defined above
for key, rows in display_shape.items():
    plot_tab = []
    data_sources[key] = []
    for cur_row_idx, cur_row in enumerate(rows):
        data_sources_rows = []
        plot_row = []
        for title, (x_axis_label, y_axis_label), plot, data_idx in cur_row:
            if plot.startswith("xy_history") or plot == "bar":
                cur_data_source = ColumnDataSource({'colors': [], 'x': [], 'y': []})
            elif plot == "donut":
                cur_data_source = ColumnDataSource({'colors': [], 'x': [], 'y': [], 'y_label': []})
            else:
                cur_data_source = ColumnDataSource({'x': [], 'y': []})
            data_sources_rows.append(cur_data_source)
            
            # Initialize formatting specifics for each plot type
            if plot.startswith('xy'):
                if plot.startswith('xy_history'):
                    y_range = (0, 100)
                else:
                    y_range = None
                cur_fig = figure(plot_width=375, plot_height=375, x_axis_type='datetime', y_range=y_range)
                tooltips = []
                x_formatter = '%m/%d/%y'
                formatters = {'@x': 'datetime'}
                cur_fig.xaxis.formatter=DatetimeTickFormatter(microseconds=["%m/%d/%y"],
                                                            milliseconds=["%m/%d/%y"],
                                                                seconds=["%m/%d/%y"],
                                                                    minsec=["%m/%d/%y"],
                                                                minutes=["%m/%d/%y"],
                                                                hourmin=["%m/%d/%y"],
                                                                    hours=["%m/%d/%y"],
                                                                    days=["%m/%d/%y"],
                                                                    months=["%m/%d/%y"],
                                                                    years=["%m/%Y"])
                tooltips += [(x_axis_label, '@x{' + x_formatter + '}'), (y_axis_label, '@y{0.000}')]
                hover = HoverTool(tooltips=tooltips, formatters=formatters)
                cur_fig.add_tools(hover)
                #Change the display angle of x axis tick labels
                cur_fig.xaxis.major_label_orientation = np.pi/4
                # Set the formatting of the data point symbol for scatter plots
                if plot.startswith('xy_history'):
                    cur_fig.circle('x', 'y', source=cur_data_source, size=12, color='colors', alpha=.8)
                else:
                    cur_fig.circle('x', 'y', source=cur_data_source, size=12, color="navy", alpha=0.5)
            elif plot == "bar":
                cur_data_source.data['x'] = ["Completeness", "False Positive Rate"]
                
                cur_fig = figure(x_range=cur_data_source.data['x'], y_range=(0, 100), plot_width=375, plot_height=375)
                cur_fig.vbar(x='x', top='y', source=cur_data_source, width=0.9, fill_color='colors', line_color="black")
                cur_fig.xaxis.major_label_text_font_size = "10pt"
                cur_fig.xgrid.grid_line_color = None
                cur_fig.y_range.start = 0
            elif plot == "donut":
                cur_fig = figure(plot_width=375, plot_height=375, title="Title", toolbar_location=None, x_range=(-1.1, 1.9), y_range=(-0.5, 2.5), outline_line_alpha=0.0)
                cur_fig.annular_wedge(x=0.4, y=0.5,  inner_radius=1.0, outer_radius=1.5, direction='anticlock',
                                    start_angle=cumsum('y', include_zero=True), end_angle=cumsum('y'),
                                    line_color="white", fill_color='colors', source=cur_data_source)
                label = LabelSet(x=0, y=0, text='y_label', level='glyph', source=cur_data_source, render_mode='canvas', text_font_size="20pt")
                cur_fig.add_layout(label)
                cur_fig.axis.axis_label=None
                cur_fig.axis.visible=False
                cur_fig.grid.grid_line_color = None
                
            cur_fig.xaxis.axis_label_text_font_size = "15px"
            cur_fig.yaxis.axis_label_text_font_size = "15px"
            cur_fig.title.text_font_size = "16px"
            cur_fig.min_border = 70

            if title:
                cur_fig.title.text = title
            if x_axis_label:
                cur_fig.xaxis.axis_label = x_axis_label
            if y_axis_label:
                cur_fig.yaxis.axis_label = y_axis_label

            plot_row.append(cur_fig)

        plot_tab.append(plot_row)
        data_sources[key].append(data_sources_rows)

    if key != "summary":
        plot_tab = [[select_bar]] + plot_tab
        child_plot = gridplot(plot_tab, plot_width=375, plot_height=375)
    else:
        summary_select_bar = Select(title="Branch:", options=['IOD', 'Detect Sidereal', 'Detect Target', 'Calibrate Sidereal'])
        summary_select_bar.on_change('value', onSummarySelectChange)

        child_plot = grid([
                summary_time_series(),
                summary_select_bar,
                leaderboard()
            ])

    tabs[key] = Panel(child=child_plot, title=makePrettyName(key))

# Build the tabs for the summary and branch-specific pages
tabs = Tabs(tabs=[ value for key, value in tabs.items()])

tabs.on_change('active', onTabChange)   

# The tabs are organized in a column layout across the top of the dashboard
layout = column(tabs)

redrawPlot(current_tab, None)
#show(layout)
curdoc().add_root(layout)
