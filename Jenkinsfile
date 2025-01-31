#!/usr/bin/env groovy
/*
 * Jenkins Pipeline for ScorpioDR, based on the DRAGONS one
 *
 */

// @Library('dragons_ci@master') _

// Change these to automatically skip steps
def runtests_scorpio = 1  // 1 to enable

pipeline {

    agent any

    //triggers {
    //   //Polls Source Code Manager every 15 mins
    //   pollSCM('*/15 * * * *')
    //}

    options {
        skipDefaultCheckout(true)
        buildDiscarder(logRotator(numToKeepStr: '5'))
        timestamps()
        timeout(time: 6, unit: 'HOURS')
    }

    environment {
        MPLBACKEND = "agg"
        PATH = "$JENKINS_CONDA_HOME/bin:$PATH"
    }

    stages {

        stage ("Prepare"){
            steps{
                echo "Step would notify STARTED when dragons_ci is available"
                // sendNotifications 'STARTED'
            }
        }


        stage('Pre-install') {
            agent { label "conda" }
            environment {
                TMPDIR = "${env.WORKSPACE}/.tmp/conda/"
            }
            steps {
                echo "Update the Conda base install for all on-line nodes"
                checkout scm
                sh '.jenkins/scripts/setup_agent.sh'
                // // This is only needed if we have parallel stages later:
                // echo "Create a trial Python 3.10 env, to cache new packages"
                // sh 'tox -e py310-noop -v -r -- --basetemp=${DRAGONS_TEST_OUT} ${TOX_ARGS}'
            }
            post {
                always {
                    echo "Deleting conda temp workspace ${env.WORKSPACE}"
                    cleanWs()
                    dir("${env.WORKSPACE}@tmp") {
                      deleteDir()
                    }
                }
            }
        }

        // // Uncomment this once we have some unit tests:
        // stage('Quicker tests') {
        //     parallel {

        //         stage('Unit tests') {

        //             agent{
        //                 label "centos7"
        //             }
        //             environment {
        //                 MPLBACKEND = "agg"
        //                 DRAGONS_TEST_OUT = "unit_tests_outputs/"
        //                 TOX_ARGS = "scorpiodr scorpio_instruments"
        //                 TMPDIR = "${env.WORKSPACE}/.tmp/unit/"
        //             }
        //             steps {
        //                 echo "Running build #${env.BUILD_ID} on ${env.NODE_NAME}"
        //                 checkout scm
        //                 sh '.jenkins/scripts/setup_dirs.sh'
        //                 echo "Running tests with Python 3.10"
        //                 sh 'tox -e py310-unit -v -r -- --basetemp=${DRAGONS_TEST_OUT} --junit-xml reports/unittests_results.xml ${TOX_ARGS}'
        //                 echo "Reportint coverage to CodeCov"
        //                 sh 'tox -e codecov -- -F unit'
        //             }
        //             post {
        //                 always {
        //                     junit (
        //                         allowEmptyResults: true,
        //                         testResults: '.tmp/py310-unit/reports/*_results.xml'
        //                     )
        //                     echo "Deleting Unit tests workspace ${env.WORKSPACE}"
        //                     cleanWs()
        //                     dir("${env.WORKSPACE}@tmp") {
        //                       deleteDir()
        //                     }
        //                 }
        // //                failure {
        // //                    echo "Archiving tests results for Unit Tests"
        // //                    sh "find ${DRAGONS_TEST_OUT} -not -name \\*.bz2 -type f -print0 | xargs -0 -n1 -P4 bzip2"
        // //                             archiveArtifacts artifacts: "${DRAGONS_TEST_OUT}/**"
        // //                }
        //             }
        //         }

        //         stage('Regression Tests') {
        //             agent { label "master" }
        //             environment {
        //                 MPLBACKEND = "agg"
        //                 DRAGONS_TEST_OUT = "regression_tests_outputs"
        //                 TOX_ARGS = "scorpiodr scorpio_instruments"
        //                 TMPDIR = "${env.WORKSPACE}/.tmp/regr/"
        //             }
        //             steps {
        //                 echo "Running build #${env.BUILD_ID} on ${env.NODE_NAME}"
        //                 checkout scm
        //                 echo "${env.PATH}"
        //                 sh '.jenkins/scripts/setup_dirs.sh'
        //                 echo "Regression tests"
        //                 sh 'tox -e py310-reg -v -r -- --basetemp=${DRAGONS_TEST_OUT} --junit-xml reports/regression_results.xml ${TOX_ARGS}'
        //                 echo "Reporting coverage"
        //                 sh 'tox -e codecov -- -F regression'
        //             } // end steps
        //             post {
        //                 always {
        //                     junit (
        //                         allowEmptyResults: true,
        //                         testResults: '.tmp/py310-reg/reports/*_results.xml'
        //                     )
        //                     echo "Deleting Regression Tests workspace ${env.WORKSPACE}"
        //                     cleanWs()
        //                     dir("${env.WORKSPACE}@tmp") {
        //                       deleteDir()
        //                     }
        //                 }
        //             } // end post
        //         }
        //     } // end parallel
        // }

        stage('Instrument tests') {
            parallel {
                stage('SCORPIO Tests') {
                    when {
                        expression { runtests_scorpio  == 1 }
                    }

                    agent { label "master" }
                    environment {
                        MPLBACKEND = "agg"
                        DRAGONS_TEST_OUT = "scorpio_tests_outputs"
                        TOX_ARGS = "scorpiodr scorpio_instruments"
                        TMPDIR = "${env.WORKSPACE}/.tmp/scorpio/"
                    }
                    steps {
                        echo "Running build #${env.BUILD_ID} on ${env.NODE_NAME}"
                        checkout scm
                        sh '.jenkins/scripts/setup_dirs.sh'
                        echo "Running tests"
                        sh 'tox -e py310-scorpio -v -- --basetemp=${DRAGONS_TEST_OUT} --junit-xml reports/scorpio_results.xml ${TOX_ARGS}'
                        echo "Reporting coverage"
                        sh 'tox -e codecov -- -F scorpio'
                    }  // end steps
                    post {
                        always {
                            echo "Running 'archivePlots' from inside SCORPIO Tests"
                            archiveArtifacts artifacts: "plots/*", allowEmptyArchive: true
                            junit (
                                allowEmptyResults: true,
                                testResults: '.tmp/py310-scorpio/reports/*_results.xml'
                            )
                            echo "Deleting SCORPIO Tests workspace ${env.WORKSPACE}"
                            cleanWs()
                            dir("${env.WORKSPACE}@tmp") {
                              deleteDir()
                            }
                        }  // end always
                    }  // end post
                }  // end stage
            } // end parallel
        }

    }
    post {
        success {
            echo "Step would notify SUCCESSFUL when dragons_ci is available"
            // sendNotifications 'SUCCESSFUL'
            echo "Step would notify SUCCESSFUL when dragons_ci is available"
            // sendNotifications 'SUCCESSFUL'
//            deleteDir() /* clean up our workspace */
        }
        failure {
            echo "Step would notify FAILED when dragons_ci is available"
            // sendNotifications 'FAILED'
            echo "Step would notify FAILED when dragons_ci is available"
            // sendNotifications 'FAILED'
//            deleteDir() /* clean up our workspace */
        }
        always {
            echo "Delete master workspace ${env.WORKSPACE}"
            cleanWs()
        }
    }
}
